//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
//   Created by Sean Jeffas, sj9ry@virginia.edu
//   Last Modified July 7, 2023
//
//
//   The purpose of this script is to check the elastic event
//   selection. It will plot the HCal data and W2 distributions
//   for H2 and He3 data for the kinematic point in question.
//
//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
#include "../../include/gen-ana.h"

double bg_low;
double bg_high;

double p_low;
double p_high;

double n_low;
double n_high;

double W2min;
double W2max;

TString kinematic;

//Function to get fits for the p/n spots on HCal
void get_np_spots(TCanvas *c, TString cfg, TH2D *hdxdy,TString config, TF1 **fit_result_x,TF1 **fit_result_y){

  //First project to x
  TH1D *hdx = hdxdy->ProjectionY();
  double x_cut = -1.0;

  //These parameters are determined by looking at the peaks on the plots by eye
  if(config == "GEN2"){
    p_low = -3.5;
    p_high = -2.2;
    
    n_low = -0.8;
    n_high = 0.5;
    
    bg_low = -4.0;
    bg_high = 3.0;

    x_cut = -1.5;
  }

  if(config == "GEN3"){
    p_low = -2.0;
    p_high = -0.8;
    
    n_low = -0.4;
    n_high = 0.4;
    
    bg_low = -4.0;
    bg_high = 3.0;

    x_cut = -1;
  }

  if(config == "GEN4" || config == "GEN4b"){
    p_low = -1.4;
    p_high = -0.9;
    
    n_low = -0.3;
    n_high = 0.1;
    
    bg_low = -4.0;
    bg_high = 1.5;

    x_cut = -0.4;
  }

  //total function = proton gaus + neutron gaus + 4th order bkgd
  double par[11];
  TF1 *p_xfunc = new TF1("p_xfunc","gaus",p_low,p_high);
  TF1 *n_xfunc = new TF1("n_xfunc","gaus",n_low,n_high);  
  TF1 *bg_xfunc = new TF1("bg_xfunc","pol4",bg_low,bg_high);
  TF1 *total_xfunc = new TF1("total_xfunc","gaus(0) + gaus(3) + pol4(6)",-4,2);
  
  //Do fit but do not plot the results
  hdx->Fit(p_xfunc,"qNR");  
  hdx->Fit(n_xfunc,"qNR+");  
  hdx->Fit(bg_xfunc,"qNR+"); 

  //Put the fit parameters into the array
  p_xfunc->GetParameters(&par[0]);
  n_xfunc->GetParameters(&par[3]);
  bg_xfunc->GetParameters(&par[6]);

  //Set the parameters in the total function using the results above
  total_xfunc->SetParameters(par);
  if(config == "GEN4") {
    total_xfunc->SetParLimits(4,-0.2,-0.1);
    total_xfunc->SetParLimits(5,0.2,0.6);
  }
  hdx->Fit(total_xfunc,"qR+"); 
  
  //Get the fit results
  total_xfunc->GetParameters(&par[0]);

  double px_mean = par[1];
  double px_sigma = par[2];
  double nx_mean = par[4];
  double nx_sigma = par[5];

  //For plotting purposes set the p/n function parameters from the total fit
  p_xfunc = new TF1("p_xfunc","gaus",-4,4);
  p_xfunc->SetParameters(&par[0]);

  n_xfunc = new TF1("n_xfunc","gaus",-4,4);
  n_xfunc->SetParameters(&par[3]);
  
  //Cut between p/n spot for the y-direction fitting
  int x_cut_bin = hdxdy->GetYaxis()->FindBin(x_cut);

  TH1D *hdy_p = hdxdy->ProjectionX("dy_p",0,x_cut_bin);
  TH1D *hdy_n = hdxdy->ProjectionX("dy_n",x_cut_bin,-1);

  //This direction is just a gaus fit
  TF1 *p_yfunc = new TF1("p_yfunc","gaus",-1,1);
  TF1 *bg_yfunc = new TF1("bg_yfunc","pol4",-1.5,1.5);
  TF1 *total_yfunc = new TF1("total_yfunc","gaus(0) + pol4(3)",-1.5,1.5);
  
  hdy_p->Fit(p_yfunc,"qNR+");  
  hdy_p->Fit(bg_yfunc,"qNR+");  
  
  p_yfunc->GetParameters(&par[0]);
  bg_yfunc->GetParameters(&par[3]);

  //Set the parameters in the total function using the results above
  total_yfunc->SetParameters(par);
  hdy_p->Fit(total_yfunc,"qR+"); 
  total_yfunc->GetParameters(&par[0]);

  double py_mean = par[1];
  double py_sigma = par[2];

  TF1 *n_yfunc = new TF1("n_yfunc","gaus",-1.5,1.5);
  hdy_n->Fit(n_yfunc,"qR+");  
  n_yfunc->GetParameters(&par[0]);
  
  double ny_mean = par[1];
  double ny_sigma = par[2];

  //Right now for He3 use 1 sigma and 2 sigma for H2
  double nsigma = 2;

  if(cfg == "H2") nsigma = 2;

  //Draw ellipses to show the fit results
  TEllipse *p_spot = new TEllipse(py_mean,px_mean,nsigma*py_sigma,nsigma*px_sigma);
  TEllipse *n_spot = new TEllipse(ny_mean,nx_mean,nsigma*ny_sigma,nsigma*nx_sigma);
  p_spot->SetFillStyle(0);
  n_spot->SetFillStyle(0);

  p_spot->SetLineWidth(4);
  n_spot->SetLineWidth(4);

  p_spot->SetLineColor(kRed);
  n_spot->SetLineColor(kRed);

  c->cd();

  gStyle->SetOptStat(0);

  hdxdy->SetTitle(kinematic + " " + cfg + " HCal Elastics;#Deltay (m);#Deltax (m)");
  hdxdy->Draw("colz");
  p_spot->Draw("same");
  if(cfg == "He3") n_spot->Draw("same");


  TPaveText *pt = new TPaveText(.55,.8,.88,.88,"ndc");
  pt->AddText("Cuts on good tracks");
  pt->AddText(Form("%g < W^{2} < %g",W2min,W2max));
  pt->SetFillColor(0);
  pt->Draw("same");

  cout<<cfg<<" data:"<<endl;
  cout<<"proton dx mean = "<<px_mean<<endl;
  cout<<"proton dx sigma = "<<px_sigma<<endl;
  cout<<"neutron dx mean = "<<nx_mean<<endl;
  cout<<"neutron dx sigma = "<<nx_sigma<<endl;

  cout<<"proton dy mean = "<<py_mean<<endl;
  cout<<"proton dy sigma = "<<py_sigma<<endl;
  cout<<"neutron dy mean = "<<ny_mean<<endl;
  cout<<"neutron dy sigma = "<<ny_sigma<<endl;
  cout<<"\n\n";
  

  //Save the fit result for use later
  *fit_result_x = total_xfunc;
  *fit_result_y = total_yfunc;
  //*fit_result = n_xfunc;
 
}



void hcal_np_fit(TString cfg = "GEN2", bool is_data = 1){

  kinematic = cfg;

  TString jmgr_He3 = "../../config/" + cfg + "_He3.cfg";
  TString jmgr_H2 = "../../config/" + cfg + "_H2.cfg";
  if(cfg == "GEN2") jmgr_H2 = "../../config/" + cfg + "_H2_SBS100.cfg";  
  if(cfg == "GEN4") jmgr_H2 = "../../config/" + cfg + "b_H2.cfg";  

  Utilities::KinConf kin_He3 = Utilities::LoadKinConfig(jmgr_He3,is_data);
  Utilities::KinConf kin_H2 = Utilities::LoadKinConfig(jmgr_H2,is_data);

  analyzed_tree *T_He3 = Utilities::LoadAnalyzedRootFiles(kin_He3,is_data,1);
  analyzed_tree *T_H2 = Utilities::LoadAnalyzedRootFiles(kin_H2,is_data,1);
  
  // elastic cut limits
  W2min = kin_He3.W2min;
  W2max = kin_He3.W2max;

  //Cuts for H2 spot
  vector<double> dx_n = kin_H2.dx_n;
  double Nsigma_dx_n = kin_H2.Nsigma_dx_n;
  vector<double> dy_n = kin_H2.dy_n;
  double Nsigma_dy_n = kin_H2.Nsigma_dy_n;
  double dxmin_H2 = dx_n[0] - dx_n[1];
  double dxmax_H2 = dx_n[0] + dx_n[1];
  double dymin_H2 = dy_n[0] - dy_n[1];
  double dymax_H2 = dy_n[0] + dy_n[1];

  double coin_min_H2 = kin_H2.coin_min;
  double coin_max_H2 = kin_H2.coin_max;

  // Cuts for He3 spot
  dx_n = kin_He3.dx_n;
  Nsigma_dx_n = kin_He3.Nsigma_dx_n;
  dy_n = kin_He3.dy_n;
  Nsigma_dy_n = kin_He3.Nsigma_dy_n;
  double dxmin_He3 = dx_n[0] - dx_n[1];
  double dxmax_He3 = dx_n[0] + dx_n[1];
  double dymin_He3 = dy_n[0] - dy_n[1];
  double dymax_He3 = dy_n[0] + dy_n[1];

  double coin_min_He3 = kin_He3.coin_min;
  double coin_max_He3 = kin_He3.coin_max;


  //Set the histograms that will be filled
  TH2D *hdxdy_nocut_H2 = new TH2D("hdxdy_nocut_H2","",150,-2,2,150,-6,6);
  TH2D *hdxdy_Wcut_H2 = new TH2D("hdxdy_Wcut_H2","",150,-2,2,150,-6,6);
  TH2D *hdxdy_coin_H2 = new TH2D("hdxdy_coin_H2","",150,-2,2,150,-6,6);
  TH1D *hW2_all_H2 = new TH1D("hW2_all_H2","",200,0,4);
  TH1D *hW2_cut_H2 = new TH1D("hW2_cut_H2","",200,0,4);

  TH2D *hdxdy_nocut_He3 = new TH2D("hdxdy_nocut_He3","",150,-2,2,150,-6,6);
  TH2D *hdxdy_Wcut_He3 = new TH2D("hdxdy_Wcut_He3","",150,-2,2,150,-6,6);
  TH2D *hdxdy_coin_He3 = new TH2D("hdxdy_coin_He3","",150,-2,2,150,-4,2.5);
  TH1D *hW2_all_He3 = new TH1D("hW2_all_He3","",200,0,4);
  TH1D *hW2_cut_He3 = new TH1D("hW2_cut_He3","",200,0,4);  
  ////////////////////////////////////////////////////////////////

  int nevent = 0;
  int maxevent = T_H2->fChain->GetEntries();

  while(nevent < maxevent){
    T_H2->GetEntry(nevent++);   

    ////// Define all the cuts we will use on the data  ////////////////
    bool good_W2 = T_H2->W2 > W2min && T_H2->W2 < W2max;
    bool good_coin_time = T_H2->coin_time > coin_min_H2 && T_H2->coin_time < coin_max_H2;
    bool good_dy_elas = T_H2->dy > dymin_H2 && T_H2->dy < dymax_H2;
    bool good_dx_elas = T_H2->dx > dxmin_H2 && T_H2->dx < dxmax_H2;
    //////////////////////////////////////////////////////////////////////

    hW2_all_H2->Fill(T_H2->W2);
    if(good_dy_elas && good_dx_elas) hW2_cut_H2->Fill(T_H2->W2);
    hdxdy_nocut_H2->Fill(T_H2->dy,T_H2->dx);
    if(good_W2) hdxdy_Wcut_H2->Fill(T_H2->dy,T_H2->dx);
    if(!is_data){
      if(good_W2) hdxdy_coin_H2->Fill(T_H2->dy,T_H2->dx);
    } 
    else {
      if(good_W2 && good_coin_time) hdxdy_coin_H2->Fill(T_H2->dy,T_H2->dx);
    }

  }

  nevent = 0;
  maxevent = T_He3->fChain->GetEntries();

  while(nevent < maxevent){
    T_He3->GetEntry(nevent++);   

    ////// Define all the cuts we will use on the data  ////////////////
    bool good_W2 = T_He3->W2 > W2min && T_He3->W2 < W2max;
    bool good_coin_time = T_He3->coin_time > coin_min_He3 && T_He3->coin_time < coin_max_He3;
    bool good_dy_elas = T_He3->dy > dymin_He3 && T_He3->dy < dymax_He3;
    bool good_dx_elas = T_He3->dx > dxmin_He3 && T_He3->dx < dxmax_He3;
    //////////////////////////////////////////////////////////////////////

    hW2_all_He3->Fill(T_He3->W2);
    if(good_dy_elas && good_dx_elas) hW2_cut_He3->Fill(T_He3->W2);
    hdxdy_nocut_He3->Fill(T_He3->dy,T_He3->dx);
    if(good_W2) hdxdy_Wcut_He3->Fill(T_He3->dy,T_He3->dx);
    if(!is_data){
      if(good_W2) hdxdy_coin_He3->Fill(T_He3->dy,T_He3->dx);
    } 
    else {
      if(good_W2 && good_coin_time) hdxdy_coin_He3->Fill(T_He3->dy,T_He3->dx);
    }

  }
  
 
  ////////////////////// ~~~~~~~Plot the results~~~~~~~~  //////////////////////

  //Get some histograms that are projections of other histograms
  TH1D *hdx_H2 = hdxdy_coin_H2->ProjectionY("");
  hdx_H2->Scale(1/hdx_H2->GetEntries());
  
  TH1D *hdx_He3 = hdxdy_coin_He3->ProjectionY("hdx_He3");
  hdx_He3->SetLineColor(kRed);
  hdx_He3->Scale(1/hdx_He3->GetEntries());

  TH1D *hdx_nocut_He3 = hdxdy_nocut_He3->ProjectionY("hdx_nocut_He3");
  TH1D *hdx_Wcut_He3 = hdxdy_Wcut_He3->ProjectionY("hdx_Wcut_He3");
  TH1D *hdx_coin_He3 = hdxdy_coin_He3->ProjectionY("hdx_coin_He3");
  TH1D *hdy_coin_He3 = hdxdy_coin_He3->ProjectionX("hdy_coin_He3",hdxdy_coin_He3->GetYaxis()->FindBin(dxmin_He3), -1);
  TH1D *hdx_coin_H2 = hdxdy_coin_H2->ProjectionY("hdx_coin_H2");
  TH1D *hdy_coin_H2 = hdxdy_coin_H2->ProjectionX("hdy_coin_H2");
  

  
  TF1 *fit_dx_H2;
  TF1 *fit_dy_H2;
  TCanvas *c1 = new TCanvas("c1","",800,1000);  
  //Analyze the spot size for H2 data and return the fit results
  get_np_spots(c1, "H2", hdxdy_coin_H2,cfg,&fit_dx_H2,&fit_dy_H2);

  TF1 *fit_dx_He3;
  TF1 *fit_dy_He3;
  TCanvas *c2 = new TCanvas("c2","",800,1000);  
  //Analyze the spot size for He3 data and return the fit results
  get_np_spots(c2, "He3", hdxdy_coin_He3,cfg,&fit_dx_He3,&fit_dy_He3);
  
  TCanvas *c3 = new TCanvas("c3","",1000,800);  
  hdx_coin_He3->Draw();
  hdx_coin_He3->SetTitle("HCal He3 Data;#Deltax;Entries");
  //hdx_Wcut_He3->Draw();
  //hdx_Wcut_He3->SetTitle("HCal He3 Data;#Deltax;Entries");
  fit_dx_He3->Draw("same");

  TCanvas *c4 = new TCanvas("c4","",1000,800);  
  hdy_coin_He3->Draw();
  hdy_coin_He3->SetTitle("HCal He3 Data;#Deltay;Entries");
  //fit_dy_He3->Draw("same");

  TCanvas *c5 = new TCanvas("c5","",1000,800);  
  hdx_coin_H2->Draw();
  hdx_coin_H2->SetTitle("HCal H2 Data;#Deltax;Entries");
  fit_dx_H2->Draw("same");

  TCanvas *c6 = new TCanvas("c6","",1000,800);  
  hdy_coin_H2->Draw();
  hdy_coin_H2->SetTitle("HCal H2 Data;#Deltay;Entries");
  fit_dy_H2->Draw("same");

}
