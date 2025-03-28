//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
//   Created by Sean Jeffas, sj9ry@virginia.edu
//   Last Modified August 7, 2024
//
//
//   The purpose of this script is to plot the W2 data for
//   simulated events.
//
//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////

#include "../../include/gen-ana.h"


void W2_sim(){

  gErrorIgnoreLevel = kError; // Ignores all ROOT warnings
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);

  TString jmgr_file_GEN2 = "../../config/GEN2_He3.cfg";
  TString jmgr_file_GEN3 = "../../config/GEN3_He3.cfg";
  TString jmgr_file_GEN4 = "../../config/GEN4_He3.cfg";
  
  Utilities::KinConf kin_info_GEN2 = Utilities::LoadKinConfig(jmgr_file_GEN2,0);
  Utilities::KinConf kin_info_GEN3 = Utilities::LoadKinConfig(jmgr_file_GEN3,0);
  Utilities::KinConf kin_info_GEN4 = Utilities::LoadKinConfig(jmgr_file_GEN4,0);

  Utilities::KinConf kin_info_data_GEN2 = Utilities::LoadKinConfig(jmgr_file_GEN2,1);
  Utilities::KinConf kin_info_data_GEN3 = Utilities::LoadKinConfig(jmgr_file_GEN3,1);
  Utilities::KinConf kin_info_data_GEN4 = Utilities::LoadKinConfig(jmgr_file_GEN4,1);
  
  analyzed_tree *T_sim_GEN2 = Utilities::LoadAnalyzedRootFiles(kin_info_GEN2,0,0);
  analyzed_tree *T_sim_GEN3 = Utilities::LoadAnalyzedRootFiles(kin_info_GEN3,0,0);
  analyzed_tree *T_sim_GEN4 = Utilities::LoadAnalyzedRootFiles(kin_info_GEN4,0,0);

  analyzed_tree *T_data_GEN2 = Utilities::LoadAnalyzedRootFiles(kin_info_data_GEN2,1,0);
  analyzed_tree *T_data_GEN3 = Utilities::LoadAnalyzedRootFiles(kin_info_data_GEN3,1,0);
  analyzed_tree *T_data_GEN4 = Utilities::LoadAnalyzedRootFiles(kin_info_data_GEN4,1,0);

  TH1F *h_GEN2 = new TH1F("h_GEN2","Simulation Invariant Mass;W^{2} (GeV^{2});Normalized Entries",100,-2,4);
  TH1F *h_GEN3 = new TH1F("h_GEN3","",100,-2,4);
  TH1F *h_GEN4 = new TH1F("h_GEN4","",100,-2,4);

  TH1F *h_data_GEN2 = new TH1F("h_data_GEN2","Data Invariant Mass;W^{2} (GeV^{2});Normalized Entries",100,-2,4);
  TH1F *h_data_GEN3 = new TH1F("h_data_GEN3","",100,-2,4);
  TH1F *h_data_GEN4 = new TH1F("h_data_GEN4","",100,-2,4);
  
  T_sim_GEN2->fChain->Draw("W2>>h_GEN2");
  T_sim_GEN3->fChain->Draw("W2>>h_GEN3");
  T_sim_GEN4->fChain->Draw("W2>>h_GEN4");

  T_data_GEN2->fChain->Draw("W2>>h_data_GEN2","nCut");
  T_data_GEN3->fChain->Draw("W2>>h_data_GEN3","nCut");
  T_data_GEN4->fChain->Draw("W2>>h_data_GEN4","nCut");
  
  h_GEN2->Scale(1.0/h_GEN2->GetMaximum());
  h_GEN3->Scale(1.0/h_GEN3->GetMaximum());
  h_GEN4->Scale(1.0/h_GEN4->GetMaximum());

  h_data_GEN2->Scale(1.0/h_data_GEN2->GetMaximum());
  h_data_GEN3->Scale(1.0/h_data_GEN3->GetMaximum());
  h_data_GEN4->Scale(1.0/h_data_GEN4->GetMaximum());

  h_GEN3->SetLineColor(kRed);
  h_GEN4->SetLineColor(kGreen);

  h_data_GEN3->SetLineColor(kRed);
  h_data_GEN4->SetLineColor(kGreen);
  
  TCanvas *c = new TCanvas("c","",800,600);
  h_GEN2->Draw("hist");
  h_GEN3->Draw("same hist");
  h_GEN4->Draw("same hist");

  
  TLegend *legend = new TLegend(.65,.35,.88,.47);
  legend->AddEntry("h_GEN2","Kin2","l");
  legend->AddEntry("h_GEN3","Kin3","l");
  legend->AddEntry("h_GEN4","Kin4","l");
  legend->SetLineColor(0);
  legend->Draw("same");

  TCanvas *c2 = new TCanvas("c2","",800,600);
  h_data_GEN2->Draw("hist");
  h_data_GEN3->Draw("same hist");
  h_data_GEN4->Draw("same hist");

  TPaveText *pt1 = new TPaveText(.64,.21,.89,.33,"ndc");
  pt1->AddText("Neutron Spot Cut");
  pt1->SetFillColor(0);
  pt1->Draw("same");
}
