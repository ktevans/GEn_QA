//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
//   Created by Sean Jeffas, sj9ry@virginia.edu
//   Last Modified July 17, 2024
//
//
//   The purpose of this script is to take analyzed data and
//   to calculate the asymmetry. It requires the output root
//   file from QuasiElastic_ana.C.
//
//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
#include <vector>
#include <iostream>

#include "../../include/gen-ana.h"

DBparse::DBInfo DBInfo;

// Load database files
void getDB(TString cfg){
  
  cout<<"Attempting to load DB File"<<endl;
  cout<<"---------------------------------------------------------------"<<endl;

   vector<DBparse::DBrequest> request = {
    {"Helicity Quality","Helicity readback good? (0/1 = bad/good)",1},
    {"Moller Quality","Moller measurements known? (0/1 = no/yes)",1}
  };

  DBInfo.cfg = cfg;
  DBInfo.var_req = request;

  DB_load(DBInfo);

  cout<<"---------------------------------------------------------------"<<endl;

}

void Asymmetry_vs_run(TString cfg)
{
  gErrorIgnoreLevel = kError; // Ignores all ROOT warnings

  // Define a clock to get macro processing time
  TStopwatch *sw = new TStopwatch(); sw->Start();
  
  // Load the kinematics and data
  TString jmgr_file = "../../config/" + cfg + "_He3.cfg";
  Utilities::KinConf kin_info = Utilities::LoadKinConfig(jmgr_file,1);
  
  // Load database
  DBInfo.W2min = kin_info.W2min;
  DBInfo.W2max = kin_info.W2max;
  DBInfo.dymax = kin_info.dymax;
  getDB(cfg);
  
  // Set the analyzed tree data
  analyzed_tree *T = Utilities::LoadAnalyzedRootFiles(kin_info,1,0);

  // Load info from the config file
  int IHWP_Flip = kin_info.IHWP_Flip;

  double W2min = kin_info.W2min;
  double W2max = kin_info.W2max;

  vector<double> dx_n = kin_info.dx_n;
  vector<double> dy_n = kin_info.dy_n;
  vector<double> dx_p = kin_info.dx_p;
  vector<double> dy_p = kin_info.dy_p;

  double n_dxmin = dx_n[0] - dx_n[1];
  double n_dxmax = dx_n[0] + dx_n[1];
  double n_dymin = dy_n[0] - dy_n[1];
  double n_dymax = dy_n[0] + dy_n[1];
  double p_dxmin = dx_p[0] - dx_p[1];
  double p_dxmax = dx_p[0] + dx_p[1];
  double p_dymin = dy_p[0] - dy_p[1];
  double p_dymax = dy_p[0] + dy_p[1];

  double coin_min = kin_info.coin_min;
  double coin_max = kin_info.coin_max;
  
  //Totals for all runs
  int Yp_n_total = 0;
  int Ym_n_total = 0;
  int ncut_n_total = 0;
  int Yp_p_total = 0;
  int Ym_p_total = 0;
  int ncut_p_total = 0;

  //Y+ and Y- neutron helicity yields
  //For individual runs
  int Yp_n = 0;
  int Ym_n = 0;
  int ncut_n = 0;
  
  //Y+ and Y- proton helicity yields
  //For individual runs
  int Yp_p = 0;
  int Ym_p = 0;
  int ncut_p = 0;

  // Setup some variables for helicity per run
  vector<double> runs;
  vector<double> A_n;
  vector<double> A_n_err;
  vector<double> A_p;
  vector<double> A_p_err;
  vector<double> p_nevents_points, n_nevents_points, A_n_nevents, A_p_nevents, A_n_err_nevents, A_p_err_nevents;

  int currentrunnum = -1;
  int plotrunnum = 10;
  int ndata_per_point = 2000;    //Used to make the plot of A over time
  int ndata_per_point_min = 1000;//Don't plot A if its less than this many events
  int ndata  = 0;
  int past_IHWP = 0;

  int nevent = 0;
  int maxevent = T->fChain->GetEntries(); // Used to loop through the tree

  //This will loop over all events and record the helicity inside the proton or
  //neutron spots. It will also record the asymmetry over a number of events to 
  //make plots of asymmetry over time
  //while(nevent < maxevent){
  while(nevent < 70000000){

    T->GetEntry(nevent++);

    ////// Define all the cuts we will use on the data  ////////////////
    bool good_hel = DBInfo.GoodHel[T->runnum] && (T->helicity == -1 || T->helicity == 1);
    bool good_moller = DBInfo.GoodMoller[T->runnum];
    bool good_W2 = T->W2 > W2min && T->W2 < W2max;
    bool good_neutron = T->dx > n_dxmin && T->dx < n_dxmax && T->dy > n_dymin && T->dy < n_dymax;
    bool good_proton = T->dx > p_dxmin && T->dx < p_dxmax && T->dy > p_dymin && T->dy < p_dymax;
    bool good_coin_time = T->coin_time > coin_min && T->coin_time < coin_max;
    //////////////////////////////////////////////////////////////////////

    int IHWP = T->IHWP;
    int helicity = T->helicity;
    helicity *= IHWP_Flip;
    
    if(!good_hel) continue;  //Remove events with bad helicity
    if(!good_W2) continue;   //Remove events outside of W2 region
    if(!good_coin_time) continue;  //Remove events outside of coincidence time

    if(past_IHWP == 0) past_IHWP = IHWP;  //initialize the IHWP state

    //If the IHWP state has changed or we have reached the desired number of events then we want to add this point to the plot
    if(IHWP != past_IHWP || ndata == ndata_per_point){      
      
      //This is used to add space between IHWP states on the plot so it looks nicer
      if(IHWP != past_IHWP) plotrunnum += 20;

      //Calculate the Asymmetries for this data subset
      double Asym_n = (Yp_n - Ym_n)*1.0/(Yp_n + Ym_n);
      double Error_n = CalcAsymErr(Yp_n, Ym_n);
      
      double Asym_p = (Yp_p - Ym_p)*1.0/(Yp_p + Ym_p);
      double Error_p = CalcAsymErr(Yp_p, Ym_p);
      
      // Add them to the array if they are above the minimum
      if(ndata > ndata_per_point_min){
	runs.push_back(plotrunnum++);
	A_n.push_back(Asym_n*100);
	A_n_err.push_back(Error_n*100);
	
	A_p.push_back(Asym_p*100);
	A_p_err.push_back(Error_p*100);

      }      

      // Set everything back to 0 for the next data point
      Yp_n = 0;
      Ym_n = 0;
      Yp_p = 0;
      Ym_p = 0;
      ncut_n = 0;
      ncut_p = 0;
      ndata = 0;      

    } //End loop over plotting data points

    if(good_neutron){ //Cut on neutron spot
      if(helicity == 1){
	Yp_n++;
	ncut_n++;
	Yp_n_total++;
	ncut_n_total++;
	ndata++;
      }
      if(helicity == -1){
	Ym_n++;
	ncut_n++;
	Ym_n_total++;
	ncut_n_total++;
	ndata++;
      }
    }
    if(good_proton){ //Cut on proton spot
      if(helicity == 1){
	Yp_p++;
	ncut_p++;
	Yp_p_total++;
	ncut_p_total++;
      }
      if(helicity == -1){
	Ym_p++;
	ncut_p++;
	Ym_p_total++;
	ncut_p_total++;
      }
    }

    //These variables are used to plot the Asymmetry vs total events in increments of 10000 events
    if(good_proton && ncut_p_total % 10000 == 0){
       
      double Asym_p = (Yp_p_total - Ym_p_total)*1.0/(Yp_p_total + Ym_p_total);
      double Error_p = CalcAsymErr(Yp_p_total, Ym_p_total);

      p_nevents_points.push_back(ncut_p_total);
      A_p_nevents.push_back(Asym_p*100);	
      A_p_err_nevents.push_back(Error_p*100);

    }

    if(good_neutron && ncut_n_total % 10000 == 0){
	
      double Asym_n = (Yp_n_total - Ym_n_total)*1.0/(Yp_n_total + Ym_n_total);
      double Error_n = CalcAsymErr(Yp_n_total, Ym_n_total);
      
      n_nevents_points.push_back(ncut_n_total);
      A_n_nevents.push_back(Asym_n*100);	
      A_n_err_nevents.push_back(Error_n*100);
	
    }

    past_IHWP = IHWP;

  }

  // Draw plots of Asymmetry vs run number
  TGraphErrors *g_A_n = new TGraphErrors(runs.size(),&runs[0],&A_n[0],0,&A_n_err[0]);
  TGraphErrors *g_A_p = new TGraphErrors(runs.size(),&runs[0],&A_p[0],0,&A_p_err[0]);

  TCanvas *c = new TCanvas("c","",800,600);
  g_A_n->Draw("AP");
  g_A_n->SetMarkerStyle(20);
  g_A_n->SetTitle("Asymmetry vs Run Number;Run Number;Asymmetry (%)");
    
  g_A_p->Draw("P same");
  g_A_p->SetMarkerStyle(21);
  g_A_p->SetMarkerColor(kRed);
  g_A_p->SetLineColor(kRed);

  //g_A_n->GetYaxis()->SetRangeUser(-15,15);
  g_A_n->GetYaxis()->SetRangeUser(-8,8);
  g_A_n->GetXaxis()->SetLimits(0,160);

  TLine *line0 = new TLine(0,0,160,0);
  line0->SetLineStyle(2);
  line0->Draw("same");
  
  TLegend *legend = new TLegend(0.6,0.79,0.89,0.89);
  legend->AddEntry(g_A_n,"(e,e'n) events","p");
  legend->AddEntry(g_A_p,"(e,e'p) events","p");
  legend->SetLineColor(0);
  legend->Draw("same");


  //Draw plots of Asymmetry vs number of events
  TGraphErrors *g_A_n_nevents = new TGraphErrors(n_nevents_points.size(),&n_nevents_points[0],&A_n_nevents[0],0,&A_n_err_nevents[0]);
  TGraphErrors *g_A_p_nevents = new TGraphErrors(p_nevents_points.size(),&p_nevents_points[0],&A_p_nevents[0],0,&A_p_err_nevents[0]);
  
  
  TCanvas *c2 = new TCanvas("c2","",800,600);
  g_A_p_nevents->Draw("AP");
  g_A_p_nevents->SetMarkerStyle(8);
  g_A_p_nevents->SetTitle("Asymmetry vs Total Events;Events Analyzed;Asymmetry (%)");
  g_A_p_nevents->SetMarkerColor(kRed);
  
    
  g_A_n_nevents->Draw("P same");
  g_A_n_nevents->SetMarkerStyle(8);

  //g_A_n_nevents->GetYaxis()->SetRangeUser(-15,15);
  //g_A_n_nevents->GetYaxis()->SetRangeUser(-8,8);
  //g_A_n_nevents->GetXaxis()->SetLimits(0,250);

  /*
  TLegend *legend = new TLegend(0.6,0.79,0.89,0.89);
  legend->AddEntry(g_A_n,"(e,e'n) events","p");
  legend->AddEntry(g_A_p,"(e,e'p) events","p");
  legend->SetLineColor(0);
  */
  legend->Draw("same");
  

  //Calculate total asymmetry and print it to the screen
  double A_n_total = (Yp_n_total - Ym_n_total)*1.0/(Yp_n_total + Ym_n_total);
  double A_p_total = (Yp_p_total - Ym_p_total)*1.0/(Yp_p_total + Ym_p_total);
  
  double Err_n = sqrt((1 - A_n_total*A_n_total)/ncut_n_total);
  double Err_p = sqrt((1 - A_p_total*A_p_total)/ncut_p_total);

  cout<<"Total Neutron Events "<<ncut_n_total<<endl;
  cout<<"Y+ "<<Yp_n_total<<endl;
  cout<<"Y- "<<Ym_n_total<<endl;
  cout<<"A = "<<A_n_total*100<<"% +/- "<<Err_n*100<<"%"<<endl;
  cout<<endl;
  cout<<"Total Proton Events "<<ncut_p_total<<endl;
  cout<<"Y+ "<<Yp_p_total<<endl;
  cout<<"Y- "<<Ym_p_total<<endl;
  cout<<"A = "<<A_p_total*100<<"% +/- "<<Err_p*100<<"%"<<endl;


}
