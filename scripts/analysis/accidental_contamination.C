//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
//   Created by Sean Jeffas, sj9ry@virginia.edu
//   Last Modified August 7, 2024
//
//
//   The purpose of this script is to use He3 data to test
//   cuts for accidental background contributions.
//
//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////

#include "../../include/gen-ana.h"

double coin_min;
double coin_max;

double coin_bg_low;
double coin_bg_hi;

double accidental_min;
double accidental_max;

double time_low;
double time_hi;

double W2min;
double W2max;

double dxmin;
double dxmax;
double dymin;
double dymax;

DBparse::DBInfo DBInfo;

// Load database files
void getDB(TString cfg){
  
  cout<<"Attempting to load DB File"<<endl;
  cout<<"---------------------------------------------------------------"<<endl;

   vector<DBparse::DBrequest> request = {
     {"Helicity Quality","Helicity readback good? (0/1 = bad/good)",1},
     {"Moller Quality","Moller measurements known? (0/1 = no/yes)",1},
     {"Asymmetry Correction","All asymmetry correction parameters",1}
  };

  DBInfo.cfg = cfg;
  DBInfo.var_req = request;

  DB_load(DBInfo);

  cout<<"---------------------------------------------------------------"<<endl;

}

void DrawLines(TVirtualPad *pad, double xmin, double xmax, Color_t color){

  pad->Update();
  pad->cd();
  TLine *l_min = new TLine(xmin,pad->GetUymin(),xmin,pad->GetUymax());
  TLine *l_max = new TLine(xmax,pad->GetUymin(),xmax,pad->GetUymax());
  l_min->SetLineColor(color);  
  l_max->SetLineColor(color);  
  l_min->SetLineWidth(2);  
  l_max->SetLineWidth(2);  
  l_min->Draw("same");    
  l_max->Draw("same");   

}


TBox *DrawBox(TVirtualPad *pad, double xmin, double xmax, Color_t color){

  pad->Update();
  pad->cd();
  TBox *box = new TBox(xmin, pad->GetUymin(), xmax, pad->GetUymax());
  box->SetFillColorAlpha(color, 0.2);
  box->Draw("same");

  return box;
}

void DrawAccidentals(TVirtualPad *pad, TH1F *hcoin, TGraphErrors *gA, TString cut_string){

  gA->SetMarkerStyle(8);
  gA->GetXaxis()->SetLabelSize(0.1);
  gA->GetYaxis()->SetLabelSize(0.1);
  hcoin->GetYaxis()->SetLabelSize(0.05);
  gA->GetXaxis()->SetTitleSize(0.10);
  gA->GetYaxis()->SetTitleSize(0.10);
  hcoin->GetYaxis()->SetTitleSize(0.05);
  gA->GetYaxis()->SetTitleOffset(0.28);
  //hcoin->GetYaxis()->SetTitleOffset(0.35);
  
  pad->cd();
  TPad *pad1 = new TPad("pad1", "pad1", 0, 0.3, 1, 1);
  TPad *pad2 = new TPad("pad2", "pad2", 0, 0, 1, 0.3);

  pad1->SetBottomMargin(0.00); // Set bottom margin for pad1
  pad2->SetTopMargin(0.00); // Set top margin for pad2
  pad2->SetBottomMargin(0.30);
  
  pad1->Draw();
  pad2->Draw();
  
  pad1->cd();
  pad1->SetGridx();
  
  hcoin->Draw();
  pad1->Update();
  
  DrawLines(pad1->cd(), coin_min, coin_max, kRed);
  DrawLines(pad1->cd(), coin_bg_low, coin_bg_hi, kBlue);
  TBox *box_red = DrawBox(pad1->cd(), accidental_min, coin_min, kRed);
  DrawBox(pad1->cd(), coin_max, accidental_max, kRed);
  TBox *box_blue = DrawBox(pad1->cd(), coin_bg_low, coin_bg_hi, kBlue);


  TLegend *legend = new TLegend(.65,.35,.88,.47);
  legend->AddEntry(box_red,"Bad Timing Cut","f");
  legend->AddEntry(box_blue,"Accidental Fraction","f");
  legend->SetLineColor(0);
  legend->Draw("same");
  
  TPaveText *pt = new TPaveText(0.65,0.47,0.88,0.70,"ndc");
  pt->AddText("Vertex & Preshower Cut");
  pt->AddText(Form("%g < W^{2} < %g",W2min,W2max));
  pt->AddText(Form("|#Deltax| < %g",dxmax));
  pt->AddText(cut_string);
  pt->SetFillColor(0);
  pt->Draw("same");
  
  pad2->cd();
  pad2->SetGridx();

  TLine *line0 = new TLine(time_low,0,time_hi,0);
  line0->SetLineStyle(2);

  gA->Draw("AP");
  line0->Draw("same");
  gA->GetXaxis()->SetLimits(time_low,time_hi);
  gA->GetYaxis()->SetRangeUser(-6,6);

}

void accidental_contamination(TString cfg = "GEN2", double W2mincut = -100, double W2maxcut = -100, double dycut = -100){

  TString kin;

  if(cfg == "GEN3")
    kin = "Kin3";
  else if(cfg == "GEN4")
    kin = "Kin4";
  else
    kin = "Kin2";
  
  gErrorIgnoreLevel = kError; // Ignores all ROOT warnings
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);

  TString jmgr_file = "../../config/" + cfg + "_He3.cfg";
  Utilities::KinConf kin_info = Utilities::LoadKinConfig(jmgr_file,1);

  vector<double> dx_n = kin_info.dx_n;
  vector<double> dy_n = kin_info.dy_n;
  dxmin = dx_n[0] - dx_n[1];
  dxmax = dx_n[0] + dx_n[1];
  dymin = dy_n[0] - dy_n[1];
  dymax = dy_n[0] + dy_n[1];
  
  // Set variables for cuts
  W2min = kin_info.W2min;
  W2max = kin_info.W2max;
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

  double dy_bg_min = kin_info.dymin;
  double dy_bg_max = kin_info.dymax;

  

  int IHWP_Flip = kin_info.IHWP_Flip;

  coin_min = kin_info.coin_min;
  coin_max = kin_info.coin_max;

  //double coin_avg = (coin_min + coin_max) / 2;
  double coin_size = coin_max - coin_min;

  coin_bg_low = coin_max + 0.7*coin_size;
  coin_bg_hi = coin_bg_low + coin_size;
  
  //coin_bg_low = coin_max + 7*kin_info.coin_time_cut[1];
  //coin_bg_hi = coin_bg_low + 2*kin_info.Nsigma_coin_time*kin_info.coin_time_cut[1];

  accidental_min = 40;
  accidental_max = 180;

  if(cfg == "GEN3"){
    accidental_min = 60;
    accidental_max = 170;
  }
  else if(cfg == "GEN4"){
    accidental_min = 60;
    accidental_max = 170;
  }

  // Set cuts limits for histogram
  time_low = 0;
  time_hi = 200;

  double coin_bin_size = coin_size;
  const int N_points = (int)((time_hi - time_low) / coin_bin_size);
  
  double coin_low[N_points], coin_hi[N_points], coin_mid[N_points];
  double A_array[N_points], A_err_array[N_points];
  double A_inel_array[N_points], A_inel_err_array[N_points];
  int N_p[N_points], N_m[N_points];
  int N_p_inel[N_points], N_m_inel[N_points];
  
  for(int i=0; i < N_points; i++){
    coin_low[i] = time_low + i*coin_bin_size;
    coin_hi[i] = coin_low[i] + coin_bin_size;
    coin_mid[i] = (coin_hi[i] + coin_low[i]) / 2;
    N_p[i] = 0;
    N_m[i] = 0;
    N_p_inel[i] = 0;
    N_m_inel[i] = 0;
  }

  analyzed_tree *T_data = Utilities::LoadAnalyzedRootFiles(kin_info,1,1);

  /////Set the histograms
  double hxmin = -4;
  double hxmax = 2.5;

  if(cfg == "GEN2"){
    hxmin = -6;
    hxmax = 3;
  }

  TH1F *hcoin_time_elas = new TH1F("hcoin_time_elas",kin + " Quasielastic Coincidence Time;Coincidence Time (ns);Entries",100,time_low,time_hi);
  TH1F *hcoin_time_inel = new TH1F("hcoin_time_inel","Inelastic Coincidence Time;Coincidence Time (ns);Entries",100,time_low,time_hi);
  TH1F *hdx_inel = new TH1F("hdx_inel","HCal Reconstructed Position;#Deltax (m);Entries",100,hxmin,hxmax);
  TH1F *hdx_inel_acc = new TH1F("hdx_inel_acc","",100,hxmin,hxmax);
  
  
  int nevent = 0;
  int maxevent = T_data->fChain->GetEntries();
  int N_p_tot = 0;
  int N_m_tot = 0;
  int N_p_tot_inel = 0;
  int N_m_tot_inel = 0;
  int N_acc = 0;
  int N_QE = 0;
  int N_acc_inel = 0;
  int N_inel = 0;

  while(nevent < maxevent){
  //while(nevent < 5000000){
    T_data->GetEntry(nevent++);   
    
    double coin_time = T_data->coin_time;

    ////// Define all the cuts we will use on the data  ////////////////
    bool good_hel = DBInfo.GoodHel[T_data->runnum] && (T_data->helicity == -1 || T_data->helicity == 1);
    bool good_moller = DBInfo.GoodMoller[T_data->runnum];
    bool good_He3 = T_data->He3Pol > 0.01;
    bool good_W2 = T_data->W2 > W2min && T_data->W2 < W2max;
    bool dy_bg_cut = T_data->dy < dy_bg_min || T_data->dy > dy_bg_max;
    bool good_dy_elas = T_data->dy > dymin && T_data->dy < dymax;
    bool good_dx_elas = T_data->dx > dxmin && T_data->dx < dxmax;
    bool good_coin_time = T_data->coin_time > coin_min && T_data->coin_time < coin_max;
    bool bad_coin_time = coin_time > coin_bg_low && coin_time < coin_bg_hi;
    //////////////////////////////////////////////////////////////////////

    if(!good_hel) continue;
    if(!good_moller || !good_He3) continue;
    if(coin_time < time_low || coin_time > time_hi) continue; //Cut outliers
    if(!good_W2) continue;

    int helicity = T_data->helicity;
    helicity *= -1*T_data->IHWP*IHWP_Flip;

    int coin_bin = (int) ((coin_time - time_low) / coin_bin_size);
    
    if(dy_bg_cut){  //This is to look at the inelastic data
      if(good_dx_elas){
	hcoin_time_inel->Fill(coin_time);
	
	if(helicity == 1){
	  N_p_inel[coin_bin]++;
	  if(!good_coin_time) N_p_tot_inel++;
	}
	else if(helicity == -1){
	  N_m_inel[coin_bin]++;
	  if(!good_coin_time) N_m_tot_inel++;
	}
      }

      if(good_coin_time){
	if(good_dx_elas) N_inel++;
	hdx_inel->Fill(T_data->dx);
      }
      if(bad_coin_time){
	if(good_dx_elas) N_acc_inel++;
	hdx_inel_acc->Fill(T_data->dx);
      }
    }

    if(!good_dx_elas) continue;
    if(!good_dy_elas) continue;

    hcoin_time_elas->Fill(coin_time);  // Fill coincidence histogram
  
    if(good_coin_time)
      N_QE++;

    if(bad_coin_time)
      N_acc++;

    if(helicity == 1){
      if(!good_coin_time)
	N_p_tot++;
      N_p[coin_bin]++;
    }
    else if(helicity == -1){
      if(!good_coin_time)
	N_m_tot++;
      N_m[coin_bin]++;
    }
  
    
  }


  // Calculate our asymmetries
  for(int ibin=0; ibin < N_points; ibin++){
    A_array[ibin] = (N_p[ibin] - N_m[ibin])*1.0/(N_p[ibin] + N_m[ibin]) * 100;
    A_err_array[ibin] = CalcAsymErr(N_p[ibin],N_m[ibin]) * 100;

    A_inel_array[ibin] = (N_p_inel[ibin] - N_m_inel[ibin])*1.0/(N_p_inel[ibin] + N_m_inel[ibin]) * 100;
    A_inel_err_array[ibin] = CalcAsymErr(N_p_inel[ibin],N_m_inel[ibin]) * 100;
  }

  double Delta = N_p_tot - N_m_tot;
  double Sigma = N_p_tot + N_m_tot;
  double A_acc = Delta / Sigma;
  double f_acc = 1.0*N_acc / N_QE;
  double Aerr = CalcAsymErr(N_p_tot,N_m_tot);
  double Ferr = CalcFractionErr(N_acc, N_QE);

  double A_acc_inel = 1.0*(N_p_tot_inel - N_m_tot_inel) / (N_p_tot_inel + N_m_tot_inel);
  double Aerr_inel = CalcAsymErr(N_p_tot_inel,N_m_tot_inel);
  double f_acc_inel = 1.0*N_acc_inel / N_inel;
  double Ferr_inel = CalcFractionErr(N_acc_inel, N_inel);

  DBInfo.AccidentalAsymmetry = A_acc;
  DBInfo.AccidentalAsymmetryErr = Aerr;
  DBInfo.AccidentalFraction = f_acc;
  DBInfo.AccidentalFractionErr = Ferr;
  
  DB_SetCorrections(DBInfo);
  
  cout<<"--------- Elastic Accidental Info ---------"<<endl;
  cout<<"Accidental Np = "<<N_p_tot<<"   and Nm = "<<N_m_tot<<endl;
  cout<<"Accidental N = "<<N_acc<<endl;
  cout<<"Accidental Asymmetry = "<<A_acc<<"  Aerr = "<<Aerr<<endl;
  cout<<"Accidental Fraction = "<<f_acc<<"  Ferr = "<<Ferr<<endl;
  

  // Make all the graphs
  TGraphErrors *gA = new TGraphErrors(N_points,coin_mid,A_array,0,A_err_array);
  TGraphErrors *gA_inel = new TGraphErrors(N_points,coin_mid,A_inel_array,0,A_inel_err_array);
  
  gA->SetTitle(";Coincidence Time (ns);A (%)");
  gA_inel->SetTitle(";Coincidence Time (ns);A (%)");

  TString cut_string_elas = Form("|#Deltay| < %g",dymax);
  TString cut_string_inel = Form("|#Deltay| > %g",dy_bg_max);

  //TCanvas *c1 = new TCanvas("c1","",1600,600);
  //c1->Divide(2,1);
  TCanvas *c1 = new TCanvas("c1","",800,600);
  DrawAccidentals(c1->cd(), hcoin_time_elas, gA, cut_string_elas);
  //DrawAccidentals(c1->cd(2), hcoin_time_inel, gA_inel, cut_string_inel);


  TCanvas *c2 = new TCanvas("c2","",800,600);
  hdx_inel_acc->SetLineColor(kRed);
  hdx_inel->Draw();
  hdx_inel_acc->Draw("same");

  TLegend *legend = new TLegend(.58,.76,.88,.88);
  legend->AddEntry(hdx_inel,"Inel Cut","l");
  legend->AddEntry(hdx_inel_acc,"Inel + Accid Cut","l");
  legend->SetLineColor(0);
  legend->Draw("same");


  TString plot_dir = "../../plots/";
  TString plot_name = "Accidental_Cont_" + cfg + ".pdf";

  TString outputfile = plot_dir + plot_name;

  c1->Print(outputfile);
  //c2->Print(outputfile + ")");

}


