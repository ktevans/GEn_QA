//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
//   Created by Sean Jeffas, sj9ry@virginia.edu
//   Last Modified July 7, 2024
//
//
//   The purpose of this script is to check the elastic event
//   selection. It will plot the HCal data and W2 distributions
//   for H2 and He3 data for the kinematic point in question.
//
//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
#include "../../include/gen-ana.h"


void hcal_spot_plot(TString cfg = "GEN2", bool is_data = 1){

  gStyle->SetOptStat(0);

  TString kin;

  if(cfg == "GEN3")
    kin = "Kin3";
  else if(cfg == "GEN4")
    kin = "Kin4";
  else
    kin = "Kin2";
  
  TString jmgr_He3 = "../../config/" + cfg + "_He3.cfg";

  Utilities::KinConf kin_He3 = Utilities::LoadKinConfig(jmgr_He3,is_data);
  analyzed_tree *T_He3 = Utilities::LoadAnalyzedRootFiles(kin_He3,is_data,1);
  
  // elastic cut limits
  double W2min = kin_He3.W2min;
  double W2max = kin_He3.W2max;

  // Cuts for He3 spot
  vector<double> dx_n = kin_He3.dx_n;
  double Nsigma_dx_n = kin_He3.Nsigma_dx_n;
  vector<double> dy_n = kin_He3.dy_n;
  double Nsigma_dy_n = kin_He3.Nsigma_dy_n;
  double dxmin_He3 = dx_n[0] - dx_n[1];
  double dxmax_He3 = dx_n[0] + dx_n[1];
  double dymin_He3 = dy_n[0] - dy_n[1];
  double dymax_He3 = dy_n[0] + dy_n[1];

  double coin_min_He3 = kin_He3.coin_min;
  double coin_max_He3 = kin_He3.coin_max;


  //Set the histograms that will be filled
  TH2D *hdxdy_He3 = new TH2D("hdxdy_He3",kin + " ^{3}He HCal Quasielastics;#Deltay (m);#Deltax (m)",150,-2,2,150,-6,6);
  
  ////////////////////////////////////////////////////////////////

  int nevent = 0;
  int maxevent = T_He3->fChain->GetEntries();

  while(nevent < maxevent){
  //while(nevent < 2000000){
    T_He3->GetEntry(nevent++);   

    ////// Define all the cuts we will use on the data  ////////////////
    bool good_W2 = T_He3->W2 > W2min && T_He3->W2 < W2max;
    bool good_coin_time = T_He3->coin_time > coin_min_He3 && T_He3->coin_time < coin_max_He3;
    bool good_dy_elas = T_He3->dy > dymin_He3 && T_He3->dy < dymax_He3;
    bool good_dx_elas = T_He3->dx > dxmin_He3 && T_He3->dx < dxmax_He3;
    //////////////////////////////////////////////////////////////////////

    if(good_W2 && good_coin_time)
      hdxdy_He3->Fill(T_He3->dy,T_He3->dx);

  }
  
 
  ////////////////////// ~~~~~~~Plot the results~~~~~~~~  //////////////////////

  TCanvas *c2 = new TCanvas("c2","",800,1000);  
  hdxdy_He3->Draw("colz");

  TPaveText *pt = new TPaveText(.501,.775,.892,.892,"ndc");
  pt->AddText("Vertex & Preshower Cut");
  pt->AddText(Form("%g < W^{2} < %g",W2min,W2max));
  pt->AddText("Coin Time Cut");
  pt->SetFillColor(0);
  pt->Draw("same");
  
  TEllipse *n_spot = new TEllipse(dy_n[0],dx_n[0],dy_n[1],dx_n[1]);
  n_spot->SetFillStyle(0);
  n_spot->SetLineWidth(4);
  n_spot->SetLineColor(kRed);
  n_spot->Draw("same");

}
