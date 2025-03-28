//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
//   Created by Sean Jeffas, sj9ry@virginia.edu
//   Last Modified August 7, 2024
//
//
//   The purpose of this script is to calculate the pion 
//   asymmetry correction.
//
//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////

#include "../../include/gen-ana.h"

TH1F *hpse_pi;
TH1F *hpse_e;

int nbins = 150;
double hxmin = 0;
double hxmax = 2;
double acc_scale = 0.0;


DBparse::DBInfo DBInfo;

// Load database files
void getDB(TString cfg){
  
  cout<<"Attempting to load DB File"<<endl;
  cout<<"---------------------------------------------------------------"<<endl;

   vector<DBparse::DBrequest> request = {
     {"Helicity Quality","Helicity readback good? (0/1 = bad/good)",1},
     {"Asymmetry Correction","All asymmetry correction parameters",1}
  };

  DBInfo.cfg = cfg;
  DBInfo.var_req = request;

  DB_load(DBInfo);

  cout<<"---------------------------------------------------------------"<<endl;

}

struct ShowerHists{

  TH1F *hdata;
  
  TH1F *hpi_final;
  TH1F *he_final;
  TH1F *htotal_final;
  TString Type;

  double Cpi;  // Fit parameter
  double Ce;  // Fit parameter
  double Cpi_err;  // Fit parameter error
  double Ce_err;  // Fit parameter error
};

// Fit for preshower simulation
double fitdist( double *x, double *par){
  double dx = x[0];
  
  double Norm_pi = par[0];
  double Norm_e = par[1];
    
  double simu = Norm_pi * hpse_pi->Interpolate(dx) + Norm_e * hpse_e->Interpolate(dx);
  
  return simu;   
}


void dofitting(ShowerHists &all_hists){

  TH1F *hdata_input = all_hists.hdata;

  double scale = 1.0/hdata_input->Integral();
  hdata_input->Scale(scale);
  
  //Set fit to function fitsim
  TF1 *FitFunc = new TF1( "FitFunc", fitdist,hxmin,hxmax,2);

  //Set some arbitrary starting values
  FitFunc->SetNpx(1000);
  double startpar[] = {1.0,1.0};
  FitFunc->SetParameters(startpar);
  FitFunc->SetParLimits(0,0.0001,10);
  FitFunc->SetParLimits(1,0.0001,10);
  
  hdata_input->Fit(FitFunc,"q0","",hxmin,hxmax);

  double R_pi = FitFunc->GetParameter(0);
  double R_e = FitFunc->GetParameter(1);
  double R_pi_err = FitFunc->GetParError(0);
  double R_e_err = FitFunc->GetParError(1);
  
  TH1F *hpse_total = new TH1F("hpse_total_" + all_hists.Type,"",nbins,hxmin,hxmax);
  TH1F *hpse_pi_fit = (TH1F*)hpse_pi->Clone("hpse_pi_" + all_hists.Type);
  TH1F *hpse_e_fit = (TH1F*)hpse_e->Clone("hpse_e_" + all_hists.Type);

  all_hists.Cpi = 1.0/scale*R_pi;
  all_hists.Ce = 1.0/scale*R_e;
  all_hists.Cpi_err = 1.0/scale*R_pi_err;
  all_hists.Ce_err = 1.0/scale*R_e_err;
  
  //xshdata_input->Scale(1.0/scale*acc_scale);
  //hpse_pi_fit->Scale(1.0/scale*R_pi*acc_scale);
  //hpse_e_fit->Scale(1.0/scale*R_e*acc_scale);

  hdata_input->Scale(1.0/scale);
  hpse_pi_fit->Scale(all_hists.Cpi);
  hpse_e_fit->Scale(all_hists.Ce);
  
  for(int ibin = 0; ibin < nbins;ibin++){
    hpse_total->SetBinContent(ibin,hpse_pi_fit->GetBinContent(ibin) + hpse_e_fit->GetBinContent(ibin));
  }
  
  all_hists.hpi_final = hpse_pi_fit;
  all_hists.he_final = hpse_e_fit;
  all_hists.htotal_final = hpse_total;
  
  all_hists.hdata->SetMarkerStyle(kFullCircle);
  all_hists.hpi_final->SetLineColor(kRed);
  all_hists.htotal_final->SetLineColor(kMagenta);

  
}


void pion_contamination(TString cfg = "GEN2", double W2mincut = -100, double W2maxcut = -100, double dycut = -100){

  TString kin;

  if(cfg == "GEN3")
    kin = "Kin3";
  else if(cfg == "GEN4")
    kin = "Kin4";
  else
    kin = "Kin2";
  
  gErrorIgnoreLevel = kError; // Ignores all ROOT warnings
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);

  

  double PS_cut = 0.2; //GeV
  
  TString jmgr_He3 = "../../config/" + cfg + "_He3.cfg";
  TString jmgr_pim = "../../config/" + cfg + "_pim_sim.cfg";
  Utilities::KinConf kin_He3 = Utilities::LoadKinConfig(jmgr_He3,1);
  Utilities::KinConf kin_pim = Utilities::LoadKinConfig(jmgr_pim,0);
  analyzed_tree *T_data = Utilities::LoadAnalyzedRootFiles(kin_He3,1,0,1);
  analyzed_tree *T_sim_He3 = Utilities::LoadAnalyzedRootFiles(kin_He3,0,0,1);
  analyzed_tree *T_sim_pim = Utilities::LoadAnalyzedRootFiles(kin_pim,0,0,1);

  int IHWP_Flip = kin_He3.IHWP_Flip;
  double W2min = kin_He3.W2min;
  double W2max = kin_He3.W2max;

  double dy_bg_min = kin_He3.dymin;
  double dy_bg_max = kin_He3.dymax;
  vector<double> dx_n = kin_He3.dx_n;
  vector<double> dy_n = kin_He3.dy_n;
  double dxmin = dx_n[0] - dx_n[1];
  double dxmax = dx_n[0] + dx_n[1];
  double dymin = dy_n[0] - dy_n[1];
  double dymax = dy_n[0] + dy_n[1];

  // Set variables for cuts
  if(W2mincut > -10 && W2maxcut > -10 && dycut > -10){
    W2min = W2mincut;
    W2max = W2maxcut;
    dymax = dycut;
    dymin = -1*dymax;
    dxmin = dymin;
    dxmax = dymax;
  }

  DBInfo.W2min = W2min;
  DBInfo.W2max = W2max;
  DBInfo.dymax = dymax;
  getDB(cfg);

  
  double coin_min = kin_He3.coin_min;
  double coin_max = kin_He3.coin_max;
  double coin_size = coin_max - coin_min;
  double coin_bg_low = coin_max + 0.7*coin_size;
  double coin_bg_hi = coin_bg_low + coin_size;
  

  TH1F *hcoin_time_low = new TH1F("hcoin_time_low","Coincidence Time;Coincidence Time (ns);Entries",100,0,200);
  TH1F *hcoin_time_hi = new TH1F("hcoin_time_hi","Coincidence Time;Coincidence Time (ns);Entries",100,0,200);
  TH1F *hcoin_time_all = new TH1F("hcoin_time_all","Coincidence Time;Coincidence Time (ns);Entries",100,0,200);
  TH1F *hpse_bad_coin = new TH1F("hpse_bad_coin","",nbins,hxmin,hxmax);

  TH1F *hpse_low = new TH1F("hpse_low","",nbins,hxmin,hxmax);
  TH1F *hpse_hi = new TH1F("hpse_hi","He3 Simulation Preshower;Preshower Energy (GeV);",nbins,hxmin,hxmax);
  TH1F *hpse_pim_low = new TH1F("hpse_pim_low","",nbins,hxmin,hxmax);
  TH1F *hpse_pim_hi = new TH1F("hpse_pim_hi","#pi^{-} Simulation Preshower;Preshower Energy (GeV);",nbins,hxmin,hxmax);


  TH1F *hpse_data = new TH1F("hpse_data",kin + " Preshower Energy;Energy (GeV);",nbins,hxmin,hxmax);
  TH1F *hpse_data_p = new TH1F("hpse_data_p",kin + " Preshower Energy, +1 Helicity;Energy (GeV);",nbins,hxmin,hxmax);
  TH1F *hpse_data_m = new TH1F("hpse_data_m",kin + " Preshower Energy, -1 Helicity;Energy (GeV);",nbins,hxmin,hxmax);

  TH1F *hpse_data_inel = new TH1F("hpse_data_inel","Preshower Energy;Energy (GeV);",nbins,hxmin,hxmax);
  TH1F *hpse_data_inel_p = new TH1F("hpse_data_p_inel","Preshower Energy, +1 Helicity;Energy (GeV);",nbins,hxmin,hxmax);
  TH1F *hpse_data_inel_m = new TH1F("hpse_data_m_inel","Preshower Energy, -1 Helicity;Energy (GeV);",nbins,hxmin,hxmax);

  hpse_pi = new TH1F("hpse_pi","",nbins,hxmin,hxmax);
  hpse_e = new TH1F("hpse_e","Preshower Energy;Energy (GeV);",nbins,hxmin,hxmax);


  int nevent = 0;
  double N_QE = 0;
  double N_acc = 0;
  double N_PS = 0;
  int maxevent = T_data->fChain->GetEntries(); // Used to loop through the tree
  
  while(nevent < maxevent){
    //while(nevent < 15000000){
    T_data->GetEntry(nevent++);  //Get data for one event

    bool good_PS = T_data->ePS>0.001;
    bool PS_QE = T_data->ePS>0.2;
    bool good_hel = DBInfo.GoodHel[T_data->runnum] && (T_data->helicity == -1 || T_data->helicity == 1);
    bool good_W2 = T_data->W2 > W2min && T_data->W2 < W2max;
    bool good_dy_elas = T_data->dy > dymin && T_data->dy < dymax;
    bool good_dx_elas = T_data->dx > dxmin && T_data->dx < dxmax;
    bool dy_bg_cut = T_data->dy < dy_bg_min || T_data->dy > dy_bg_max;
    bool good_coin_time = T_data->coin_time > coin_min && T_data->coin_time < coin_max;
    bool bad_coin_time = T_data->coin_time > coin_bg_low && T_data->coin_time < coin_bg_hi;

    int helicity = T_data->helicity;
    helicity *= -1*T_data->IHWP*IHWP_Flip; 

    if(!good_PS || !good_hel) continue;

    if(good_coin_time && good_dx_elas && dy_bg_cut){
      if(helicity == 1)      
	hpse_data_inel_p->Fill(T_data->ePS);    
      if(helicity == -1)      
	hpse_data_inel_m->Fill(T_data->ePS);
    }
    
    if(!good_W2) continue;

    if(helicity == 1)
      hpse_data_p->Fill(T_data->ePS);
    if(helicity == -1)
      hpse_data_m->Fill(T_data->ePS);
    
    if(good_coin_time && good_dx_elas && dy_bg_cut)
      hpse_data_inel->Fill(T_data->ePS);       
    
    if(good_dx_elas && good_dy_elas){
      if(good_coin_time) hpse_data->Fill(T_data->ePS);
      if(bad_coin_time) hpse_bad_coin->Fill(T_data->ePS);
      if(T_data->ePS < PS_cut) hcoin_time_low->Fill(T_data->coin_time);
      if(T_data->ePS > PS_cut) hcoin_time_hi->Fill(T_data->coin_time);
      hcoin_time_all->Fill(T_data->coin_time);
      if(good_coin_time) N_PS += 1;
      if(bad_coin_time) N_acc += 1;
      if(PS_QE) N_QE += 1;
    }
    
  }

  TCut simCut = Form("(ePS > 0.01 && W2 > %g && W2 < %g)*weight",W2min,W2max);
  TCut hiPSCut = "(ePS > 0.01 && trP > 2.7)*weight";
  TCut lowPSCut = "(ePS > 0.01 && trP < 2.7)*weight";

  // Plots for PS dist with differrent momentum cuts
  T_sim_He3->fChain->Draw("ePS>>hpse_low",lowPSCut);
  T_sim_He3->fChain->Draw("ePS>>hpse_hi",hiPSCut);
  T_sim_pim->fChain->Draw("ePS>>hpse_pim_low",lowPSCut);
  T_sim_pim->fChain->Draw("ePS>>hpse_pim_hi",hiPSCut);

  T_sim_pim->fChain->Draw("ePS>>hpse_pi",simCut);
  T_sim_He3->fChain->Draw("ePS>>hpse_e",simCut);
  
  double scale_pi = 1.0/hpse_pi->Integral();
  double scale_e = 1.0/hpse_e->Integral();

  hpse_pi->Scale(scale_pi);
  hpse_e->Scale(scale_e);
  hpse_hi->Scale(1.0 / hpse_hi->Integral());
  hpse_low->Scale(1.0 / hpse_low->Integral());
  hpse_pim_hi->Scale(1.0 / hpse_pim_hi->Integral());
  hpse_pim_low->Scale(1.0 / hpse_pim_low->Integral());

  double S_pi_all = GetYield(hpse_pi,0,hxmax);
  double S_pi_cut = GetYield(hpse_pi,PS_cut,hxmax);
  
  
  double f_acc = 1.0*N_acc / N_PS;
  double f_acc_err = CalcFractionErr(N_acc, N_PS);
  //acc_scale = 1.0 - f_acc;

  TH1F *hsub = (TH1F*)hpse_data->Clone("hsub");
  hsub->Add(hpse_bad_coin,-1);  

  ShowerHists fit_p, fit_m, fit_all, fit_inel_p, fit_inel_m, fit_all_inel;
  fit_p.hdata = hpse_data_p;
  fit_p.Type = "p";
  fit_m.hdata = hpse_data_m;
  fit_m.Type = "m";
  fit_all.hdata = hpse_data;
  fit_all.Type = "all";
  fit_all_inel.hdata = hpse_data_inel;
  fit_all_inel.Type = "all_inel";
  fit_inel_p.hdata = hpse_data_inel_p;
  fit_inel_p.Type = "p_inel";
  fit_inel_m.hdata = hpse_data_inel_m;
  fit_inel_m.Type = "m_inel";

  dofitting(fit_p);
  dofitting(fit_m);
  dofitting(fit_all);
  dofitting(fit_all_inel);
  dofitting(fit_inel_p);
  dofitting(fit_inel_m);

  TCanvas *c1 = new TCanvas("c1","",1600,600);
  c1->Divide(2,1);

  c1->cd(1);
  hpse_hi->SetLineColor(kRed);
  hpse_hi->Draw("hist");
  hpse_low->Draw("same hist");

  TLegend *legend = new TLegend(0.65,0.72,0.89,0.89);
  legend->AddEntry("hpse_hi","Track p > 2.7 GeV","l");
  legend->AddEntry("hpse_low","Track p < 2.7 GeV","l");
  legend->SetLineColor(0);
  legend->Draw("same");

  c1->cd(2);
  hpse_pim_hi->SetLineColor(kRed);
  hpse_pim_hi->Draw("hist");
  hpse_pim_low->Draw("same hist");

  
  TCanvas *c2 = new TCanvas("c2","",1600,600);
  c2->Divide(2,1);
  
  c2->cd(1);
  fit_p.hdata->Draw();
  fit_p.hpi_final->Draw("same hist");
  fit_p.he_final->Draw("same hist");
  fit_p.htotal_final->Draw("same hist");
  
  TLegend *legend2 = new TLegend(0.59,0.65,0.89,0.89);
  legend2->AddEntry(fit_p.hdata->GetName(),"He3 Data","p");
  legend2->AddEntry(fit_p.htotal_final->GetName(),"Total Fit","l");
  legend2->AddEntry(fit_p.hpi_final->GetName(),"#pi^{-} Simulation","l");
  legend2->AddEntry(fit_p.he_final->GetName(),"e^{-} Simulation","l");
  legend2->SetLineColor(0);
  legend2->Draw("same");

  TPaveText *pt = new TPaveText(.65,.55,.88,.64,"ndc");
  pt->AddText("Vertex Cut");
  pt->AddText(Form("%g < W^{2} < %g",W2min,W2max)); 
  pt->SetFillColor(0);
  pt->Draw("same");

  c2->cd(2);
  fit_m.hdata->Draw();
  fit_m.hpi_final->Draw("same hist");
  fit_m.he_final->Draw("same hist");
  fit_m.htotal_final->Draw("same hist");
    
  TCanvas *c6 = new TCanvas("c6","",1600,600);
  c6->Divide(2,1);
  //hcoin_time_all->Scale(1.0/hcoin_time_all->Integral());
  //hcoin_time_low->Scale(1.0/hcoin_time_low->Integral());
  //hcoin_time_hi->Scale(1.0/hcoin_time_hi->Integral());

  c6->cd(1);
  hcoin_time_low->SetLineColor(kRed);
  hcoin_time_hi->SetLineColor(kGreen);
  hcoin_time_all->Draw("hist");
  //hcoin_time_all->Draw("hist same");
  //hcoin_time_low->Draw("hist same");
  Utilities::DrawLines(c6->cd(1), coin_bg_low, coin_bg_hi, kBlue);
  TBox *box_blue = Utilities::DrawBox(c6->cd(1), coin_bg_low, coin_bg_hi, kBlue);

  TLegend *legend61 = new TLegend(0.550,0.773,0.89,0.89);
  legend61->AddEntry(box_blue,"Accidental Fraction","f");
  legend61->SetLineColor(0);
  legend61->Draw("same");

  c6->cd(2);
  hpse_data->Draw("hist");
  hpse_bad_coin->SetLineColor(kRed);
  hpse_bad_coin->Draw("same");

  TLegend *legend62 = new TLegend(0.564,0.774,0.89,0.89);
  legend62->AddEntry(hpse_data,"Good Coin Cut","l");
  legend62->AddEntry(hpse_bad_coin,"Bad Coin Cut","l");
  legend62->SetLineColor(0);
  legend62->Draw("same");

  TCanvas *c3 = new TCanvas("c3","",800,600);
  fit_all.hdata->Draw();
  fit_all.hpi_final->Draw("same hist");
  fit_all.he_final->Draw("same hist");
  fit_all.htotal_final->Draw("same hist");

  legend2->Draw("same");

  TPaveText *pt2 = new TPaveText(.65,.55,.88,.64,"ndc");
  pt2->AddText("All QE Cuts");
  pt2->SetFillColor(0);
  pt2->Draw("same");
  

  TCanvas *c4 = new TCanvas("c4","",1600,600);
  c4->Divide(2,1);

  c4->cd(1);
  fit_inel_p.hdata->Draw();
  fit_inel_p.hpi_final->Draw("same hist");
  fit_inel_p.he_final->Draw("same hist");
  fit_inel_p.htotal_final->Draw("same hist");

  c4->cd(2);
  fit_inel_m.hdata->Draw();
  fit_inel_m.hpi_final->Draw("same hist");
  fit_inel_m.he_final->Draw("same hist");
  fit_inel_m.htotal_final->Draw("same hist");
  
  TCanvas *c5 = new TCanvas("c5","",800,600);
  fit_all_inel.hdata->Draw();
  fit_all_inel.hpi_final->Draw("same hist");
  fit_all_inel.he_final->Draw("same hist");
  fit_all_inel.htotal_final->Draw("same hist");
  
  double N_pi_p = fit_p.Cpi*S_pi_all;
  double N_pi_m = fit_m.Cpi*S_pi_all;

  double N_data_all = GetYield(fit_all.hdata,PS_cut,hxmax);
  double N_pi_all = fit_all.Cpi*S_pi_cut;

  double A = (fit_p.Cpi - fit_m.Cpi) / (fit_p.Cpi + fit_m.Cpi);
  double F = fit_all.Cpi*S_pi_cut / N_QE * (1 - 1.0*N_acc / N_PS);
    
  double Aerr = sqrt( 4.0 / pow(fit_p.Cpi + fit_m.Cpi,4) * ( pow(fit_p.Cpi*fit_m.Cpi_err,2) + pow(fit_m.Cpi*fit_p.Cpi_err,2) ));
  double Ferr = sqrt(pow(F / fit_all.Cpi * fit_all.Cpi_err,2) + pow(fit_all.Cpi*S_pi_cut /(N_QE*N_PS ),2) * N_acc + pow(fit_all.Cpi*S_pi_cut*N_acc /(N_QE*N_PS*N_PS ),2) * N_PS + F*F / N_QE);

  DBInfo.PionAsymmetry = A;
  DBInfo.PionAsymmetryErr = Aerr;
  DBInfo.PionFraction = F;
  DBInfo.PionFractionErr = Ferr;
  
  DB_SetCorrections(DBInfo);
  
  cout<<"--------------- Pion Info ---------------"<<endl;
  cout<<"Spi_all = "<<S_pi_all<<endl;
  cout<<"Spi_PS = "<<S_pi_cut<<endl;
  cout<<"Cp = "<<fit_p.Cpi<<"  Cp err = "<<fit_p.Cpi_err<<endl;
  cout<<"Cm = "<<fit_m.Cpi<<"  Cm err = "<<fit_m.Cpi_err<<endl;
  cout<<"Call = "<<fit_all.Cpi<<"  Call err = "<<fit_all.Cpi_err<<endl;
  cout<<"-----------------------------------------"<<endl;
  cout<<"Np = "<<N_pi_p<<endl;
  cout<<"Nm = "<<N_pi_m<<endl;
  cout<<"N acc = "<<N_acc<<endl;
  cout<<"N PS = "<<N_PS<<endl;
  cout<<"N pion = "<<N_pi_all<<endl;
  cout<<"-----------------------------------------"<<endl;
  cout<<"Pion Asymmetry = "<<A<<"  Aerr = "<<Aerr<<endl;
  cout<<"Pion Fraciton = "<<F<<"  Ferr = "<<Ferr<<endl;
 
      
  
  TString plot_dir = "../../plots/";
  TString plot_name = "Pion_cont_" + cfg + ".pdf";

  TString outputfile = plot_dir + plot_name;

  c1->Print(outputfile + "(");
  c6->Print(outputfile);        
  c2->Print(outputfile);
  c3->Print(outputfile + ")");
  //c4->Print(outputfile);
  //c5->Print(outputfile + ")");    
   
}


