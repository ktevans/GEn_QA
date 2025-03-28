//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
//   Created by Sean Jeffas, sj9ry@virginia.edu
//   Last Modified August 8, 2024
//
//
//   The purpose of this script is to compare real data and
//   simulated data for the same kinematic point
//
//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////

#include "../../include/gen-ana.h"


DBparse::DBInfo DBInfo;

// Load database files
void getDB(TString cfg){
  
  cout<<"Attempting to load DB File"<<endl;
  cout<<"---------------------------------------------------------------"<<endl;

   vector<DBparse::DBrequest> request = {
    {"Beam Polarization","Beam Polarization values",1},
    {"Helicity Quality","Helicity readback good? (0/1 = bad/good)",1},
    {"Moller Quality","Moller measurements known? (0/1 = no/yes)",1},
    {"Asymmetry Correction","All asymmetry correction parameters",1}
  };

  DBInfo.cfg = cfg;
  DBInfo.var_req = request;

  DB_load(DBInfo);

  cout<<"---------------------------------------------------------------"<<endl;

}

void Data_sim_compare(TString cfg = "GEN2"){

  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);

  getDB(cfg);

  distribution_fits *dists = new distribution_fits();

  if(cfg == "GEN2") dists->SetBgShapeOption("pol2");
  else dists->SetBgShapeOption("from data");
  dists->SetBgShapeOption("from data");

  TString jmgr_file = "../../config/" + cfg + "_He3.cfg";

  Utilities::KinConf kin_info = Utilities::LoadKinConfig(jmgr_file,1);

  analyzed_tree *T_data = Utilities::LoadAnalyzedRootFiles(kin_info,1,0);
  analyzed_tree *T_sim = Utilities::LoadAnalyzedRootFiles(kin_info,0,0);

  // elastic cut limits
  double W2min = kin_info.W2min;
  double W2max = kin_info.W2max;
    
  vector<double> dx_n = kin_info.dx_n;
  double Nsigma_dx_n = kin_info.Nsigma_dx_n;
  vector<double> dy_n = kin_info.dy_n;
  double Nsigma_dy_n = kin_info.Nsigma_dy_n;
  double dxmin = dx_n[0] - dx_n[1];
  double dxmax = dx_n[0] + dx_n[1];
  double dymin = dy_n[0] - dy_n[1];
  double dymax = dy_n[0] + dy_n[1];

  double dy_bg_min = kin_info.dymin;
  double dy_bg_max = kin_info.dymax;

  double coin_min = kin_info.coin_time_cut[0] - kin_info.Nsigma_coin_time*kin_info.coin_time_cut[1];
  double coin_max = kin_info.coin_time_cut[0] + kin_info.Nsigma_coin_time*kin_info.coin_time_cut[1];
  
  
  /////Set the histograms
  int nbins = 100;
  double xmin = -4;
  double xmax = 2.5;

  if(cfg == "GEN2"){
    xmin = -6;
    xmax = 3;
  }
  
  //dx
  TH1F *hdx_data = new TH1F("hdx_data","",nbins,xmin,xmax);
  TH1F *hdx_sim_p = new TH1F("hdx_sim_p","",nbins,xmin,xmax);
  TH1F *hdx_sim_n = new TH1F("hdx_sim_n","",nbins,xmin,xmax);
  TH1F *hdx_bg_data = new TH1F("hdx_bg_data","",nbins,xmin,xmax);

  
  TCut CutSimP = Form("(W2 > %g && W2 < %g && dy > %g && dy < %g && fnucl == 1) * weight",W2min,W2max,dymin,dymax);
  TCut CutSimN = Form("(W2 > %g && W2 < %g && dy > %g && dy < %g && fnucl == 0) * weight",W2min,W2max,dymin,dymax);
  
  T_sim->fChain->Draw("dx>>hdx_sim_p",CutSimP);
  T_sim->fChain->Draw("dx>>hdx_sim_n",CutSimN);

  int nevent = 0;
  int maxevent = T_data->fChain->GetEntries();

  while(nevent < maxevent){
    T_data->GetEntry(nevent++);

    ////// Define all the cuts we will use on the data  ////////////////
    bool good_hel = DBInfo.GoodHel[T_data->runnum] && (T_data->helicity == -1 || T_data->helicity == 1);
    bool good_moller = DBInfo.GoodMoller[T_data->runnum];
    bool good_He3 = T_data->He3Pol > 0.01;
    bool good_W2 = T_data->W2 > W2min && T_data->W2 < W2max;
    bool dy_bg_cut = T_data->dy < dy_bg_min || T_data->dy > dy_bg_max;
    bool good_dy_elas = T_data->dy > dymin && T_data->dy < dymax;
    bool good_dx_elas = T_data->dx > dxmin && T_data->dx < dxmax;
    bool good_coin_time = T_data->coin_time > coin_min && T_data->coin_time < coin_max;
    //////////////////////////////////////////////////////////////////////+);  

    if(!good_hel) continue;  //Remove events with bad helicity
    if(!good_moller || !good_He3) continue;  //Remove runs with a bad moller measurement or bad He3 measurements
    if(!good_coin_time) continue; 

    if(dy_bg_cut) hdx_bg_data->Fill(T_data->dx);
    if(good_W2 && good_dy_elas) hdx_data->Fill(T_data->dx);
    
  }
  
  dists->SetDataShape(hdx_data);
  dists->SetPShape(hdx_sim_p);
  dists->SetNShape(hdx_sim_n);
  dists->SetBgShape(hdx_bg_data);
  
  dists->He3_fit_dists();
  
  //Copy all the result histograms
  TH1F *hdx_data_plot = dists->GetDataHist();
  TH1F *hdx_sim_p_plot = dists->GetPHist();
  TH1F *hdx_sim_n_plot = dists->GetNHist();
  TH1F *hdx_bg_plot = dists->GetBgHist();
  TH1F *hdx_total_fit_plot = dists->GetTotalHist();

  
  gStyle->SetOptFit(0);
  
  hdx_data_plot->SetTitle("Data/Simulation Comparisons;#Deltax (m);Entries");

  hdx_data_plot->SetMarkerStyle(kFullCircle);
  hdx_total_fit_plot->SetFillColorAlpha(30,0.5);
  hdx_sim_p_plot->SetFillColorAlpha(kRed,0.3);
  hdx_sim_n_plot->SetFillColorAlpha(kBlue,0.3);
  hdx_bg_plot->SetFillColorAlpha(kMagenta,0.3);

  hdx_total_fit_plot->SetLineStyle(7);
  hdx_sim_p_plot->SetLineStyle(7);
  hdx_sim_n_plot->SetLineStyle(7);
  hdx_bg_plot->SetLineStyle(7);
  
  hdx_total_fit_plot->SetLineColor(30);
  hdx_sim_p_plot->SetLineColor(kRed);
  hdx_sim_n_plot->SetLineColor(kBlue);
  hdx_bg_plot->SetLineColor(kMagenta);
  
  TCanvas *c = new TCanvas("c","",800,600);
  hdx_data_plot->Draw();
  hdx_total_fit_plot->Draw("same hist");
  hdx_sim_p_plot->Draw("same hist");
  hdx_sim_n_plot->Draw("same hist");
  hdx_bg_plot->Draw("same hist");


  TLegend *legend = new TLegend(0.65,0.72,0.89,0.89);
  legend->AddEntry("hdx_data","Data","p");
  legend->AddEntry("hdx_total_fit","MC Fit","lf");
  legend->AddEntry("hdx_sim_p","MC p","lf");
  legend->AddEntry("hdx_sim_n","MC n","lf");
  legend->AddEntry("hdx_bg","Background","lf");
  legend->SetLineColor(0);
  legend->Draw("same");

  TPaveText *pt = new TPaveText(.65,.50,.88,.70,"ndc");
  pt->AddText("Data Background")->SetTextColor(kMagenta);
  pt->AddText("All QE Cuts");
  /*
  pt->AddText("Cuts on good tracks");
  pt->AddText("Coincidence Cuts");
  pt->AddText(Form("%g < W^{2} < %g",W2min,W2max));
  if(use_dy_cut) pt->AddText(Form("%g < #Deltay < %g",dymin,dymax));
  pt->AddText(Form("%i Neutrons",(int)hdx_sim_n_plot->GetSumOfWeights()))->SetTextColor(kBlue);
  */
  pt->SetFillColor(0);
  pt->Draw("same");


  TString output = "Data_sim_"+cfg+".pdf";
  
  c->SaveAs("../../plots/" + output);
  
}
