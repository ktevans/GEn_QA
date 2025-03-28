#include "../include/dist.h"



double distribution_fits::fitbg_pol4( double *x, double *par ){
  double dx = x[0];

  // Let's use (up to) 4th order polynomial for the background:
  double bg = 0.0;
  for( int i = 0; i<5; i++ ){
    bg += par[i]*pow(dx,i);
  }

  return bg; 
}

double distribution_fits::fitbg_pol3( double *x, double *par ){
  double dx = x[0];

  // Let's use (up to) 3rd order polynomial for the background:
  double bg = 0.0;
  for( int i = 0; i<4; i++ ){
    bg += par[i]*pow(dx,i);
  }

  return bg; 
}

double distribution_fits::fitbg_gaus( double *x, double *par ){
  double dx = x[0];

  // Let's use gaussian for the background:
  double bg = par[0]*exp(-0.5*pow((dx-par[1])/par[2],2));

  return bg; 
}

double distribution_fits::fitsim( double *x, double *par){
  double dx = x[0];

  double Norm_overall = par[0];
  double R_pn = par[1];
  double Bg_norm = par[2];

  double bg = 0.0;
  
  // Get the background function:
  if(bg_shape_option == "pol4") bg = fitbg_pol4(x,&par[3]);
  else if(bg_shape_option == "pol3") bg = fitbg_pol3(x,&par[3]);
  else if(bg_shape_option == "gaus") bg = fitbg_gaus(x,&par[3]);
  else if(bg_shape_option == "from data") bg = hdx_bg_data->Interpolate(dx);
  
  
  double simu = Norm_overall * (hdx_sim_p->Interpolate(dx) + R_pn * hdx_sim_n->Interpolate(dx) + Bg_norm * bg);
  return simu;   
}




void distribution_fits::He3_sim_fit(){

  ////////////First do some checks to make sure that things are set///////
  ////////////////////////////////////////////////////////////////////////
  // Check bg shape is set
  if(!bg_shape_set){
    cout<<"Error: [distribution_fits::He3_sim_fit] bg shape has not been set!!!"<<endl;
    exit(0);
  }
  else {
    if(bg_shape_option != "pol3" && bg_shape_option != "pol4" && bg_shape_option != "gaus" && bg_shape_option != "from data"){
      cout<<"Error: [distribution_fits::He3_sim_fit] bg shape option " + bg_shape_option + " is not supported!!!"<<endl;
      exit(0);
    }
  }

  //Check if all the histograms are set
  if(hdx_data == NULL){
    cout<<"Error: [distribution_fits::He3_sim_fit] hdx_data has not been set!!!"<<endl;
    exit(0);
  }
  if(hdx_sim_p == NULL){
    cout<<"Error: [distribution_fits::He3_sim_fit] hdx_sim_p has not been set!!!"<<endl;
    exit(0);
  }
  if(hdx_sim_n == NULL){
    cout<<"Error: [distribution_fits::He3_sim_fit] hdx_sim_n has not been set!!!"<<endl;
    exit(0);
  }
  if(hdx_bg_data == NULL){
    cout<<"Error: [distribution_fits::He3_sim_fit] hdx_bg_data has not been set!!!"<<endl;
    exit(0);
  }

  ////////////////////////////////////////////////////////////////////////
  
  double scale = hdx_data->Integral();
  hdx_data->Scale(1.0/hdx_data->Integral());
  hdx_sim_p->Scale(1.0/hdx_sim_p->Integral());
  hdx_sim_n->Scale(1.0/hdx_sim_n->Integral());
  hdx_bg_data->Scale(1.0/hdx_bg_data->Integral());

  //Get histogram data so all future histograms match
  int nbins = hdx_data->GetNbinsX();
  double xmin = hdx_data->GetXaxis()->GetBinLowEdge(1);
  double xmax = hdx_data->GetXaxis()->GetBinUpEdge(nbins);

  //Set background type and check if the option exists
  int npar = -1;
  if(bg_shape_option == "pol4") npar = 3 + 5;
  else if(bg_shape_option == "pol3") npar = 3 + 4;
  else if(bg_shape_option == "gaus") npar = 3 + 3;
  else if(bg_shape_option == "from data") npar = 3;  //Our only bg par is the scale of the bg

  //Set fit to function fitsim
  TF1 *FitFunc = new TF1( "FitFunc", this, &distribution_fits::fitsim,0,10,3);
      
  //Set some arbitrary starting values, should not be hardcoding this
  FitFunc->SetNpx(1000);
  double startpar[] = {1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0};
  FitFunc->SetParameters(startpar);
  FitFunc->SetParLimits(0,0.1,100);
  FitFunc->SetParLimits(1,0.1,100);
  FitFunc->SetParLimits(2,0,100);
    
  hdx_data->Fit(FitFunc,"0","",xmin,xmax);
    
  //Scale the p/n functions
  //This assumes our fit was like N*(p_dist + R*n_dist + N_bg * bg_dist)
  hdx_sim_p->Scale( FitFunc->GetParameter(0) );
  hdx_sim_n->Scale( FitFunc->GetParameter(0)*FitFunc->GetParameter(1) );

  hdx_bg_fit = new TH1F("hdx_bg","",nbins,xmin,xmax);
  hdx_total_fit = new TH1F("hdx_total_fit","",nbins,xmin,xmax);

  //Set the backgrond function using the fit we already found
  TF1 *fit_bg;

  if(bg_shape_option == "pol4") 
    fit_bg = new TF1("fit_bg",this, &distribution_fits::fitbg_pol4,xmin,xmax,npar - 3);
  else if(bg_shape_option == "gaus") 
    fit_bg = new TF1("fit_bg",this, &distribution_fits::fitbg_gaus,xmin,xmax,npar - 3);
  else
    fit_bg = new TF1("fit_bg",this, &distribution_fits::fitbg_pol3,xmin,xmax,npar - 3);

  for(int i=0; i<npar - 3; i++)
    fit_bg->SetParameter(i,FitFunc->GetParameter(i+3));
    
  //Fill the bg hist and total hist
  for(int ibin = 0; ibin < nbins;ibin++){
    //Again normalize bg assuming our fit is N*(p_dist + R*n_dist + N_bg * bg_dist)

    if(bg_shape_option == "from data"){
      hdx_bg_fit->SetBinContent(ibin,FitFunc->GetParameter(0) * FitFunc->GetParameter(2) * hdx_bg_data->GetBinContent(ibin));
    }
    else {
      hdx_bg_fit->SetBinContent(ibin,FitFunc->GetParameter(0) * FitFunc->GetParameter(2) * fit_bg->Eval(hdx_total_fit->GetBinCenter(ibin)));
    }
      
	
    hdx_total_fit->SetBinContent(ibin,hdx_sim_n->GetBinContent(ibin) + hdx_sim_p->GetBinContent(ibin) + hdx_bg_fit->GetBinContent(ibin));
  }


  hdx_data->Scale(scale);
  hdx_sim_p->Scale(scale);
  hdx_sim_n->Scale(scale);
  hdx_bg_fit->Scale(scale);
  hdx_total_fit->Scale(scale);

}
