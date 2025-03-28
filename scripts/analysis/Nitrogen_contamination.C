//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
//   Created by Sean Jeffas, sj9ry@virginia.edu
//   Last Modified August 7, 2024
//
//
//   The purpose of this script is to calculate the nitrogen
//   correction.
//
//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////

#include "../../include/gen-ana.h"

const int nfoil = 7;
const int nxsieve = 13;
const int nysieve = 7;
double foilz[nfoil] = {-0.225,-0.15,-0.075,0.0,0.075,0.15,0.225};

vector <Double_t> xs_cent{-(0.3+0.0492)+0.0493/cos(18.*3.14/180.),
      -(0.3+0.0492)+(0.0493+0.0492)/cos(18.*3.14/180),
      -(0.3+0.0492)+0.1493/cos(9.*3.14/180.),
      -(0.3+0.0492)+(0.1493+0.0492)/cos(9.*3.14/180.),
      -(0.3+0.0492)+(0.1493+0.0492*2.)/cos(9.*3.14/180.),
      -0.0492,
      0.0,
      0.0492,
      0.3+0.0492-(0.1493+0.0492*2.)/cos(9.*3.14/180.),
      0.3+0.0492-(0.1493+0.0492)/cos(9.*3.14/180.),
      0.3+0.0492-0.1493/cos(9.*3.14/180.),
      0.3+0.0492-(0.0493+0.0492)/cos(18.*3.14/180),
      0.3+0.0492-0.0493/cos(18.*3.14/180.)};

vector <Double_t> ys_cent;

double foil_cut = 0.021; //m


DBparse::DBInfo DBInfo;

// Load database files
void getDB(TString cfg){
  
  cout<<"Attempting to load DB File"<<endl;
  cout<<"---------------------------------------------------------------"<<endl;

   vector<DBparse::DBrequest> request = {
     {"Asymmetry Correction","All asymmetry correction parameters",1}
  };

  DBInfo.cfg = cfg;
  DBInfo.var_req = request;

  DB_load(DBInfo);

  cout<<"---------------------------------------------------------------"<<endl;

}

void DrawLines(TVirtualPad *pad, double xmin, double xmax){

  pad->Update();
  pad->cd();
  TLine *l_min = new TLine(xmin,pad->GetUymin(),xmin,pad->GetUymax());
  TLine *l_max = new TLine(xmax,pad->GetUymin(),xmax,pad->GetUymax());
  l_min->SetLineColor(kRed);  
  l_max->SetLineColor(kRed);  
  l_min->SetLineWidth(3);  
  l_max->SetLineWidth(3);  
  l_min->Draw("same");    
  l_max->Draw("same");   

}

void Nitrogen_contamination(TString cfg = "GEN2", double W2mincut = -100, double W2maxcut = -100, double dycut = -100){

  gStyle->SetOptStat(0);
  
  // y sieve must be initialized
  for (Int_t nys=0;nys<nysieve;nys++) {
    Double_t pos=nys*0.0381-0.0381*3;//old sieve
    ys_cent.push_back(pos);
  }

  double z0 = 1.25539802;
  double charge_He3 = 48.517;
  double charge_C = 0.4;

  if(cfg == "GEN4"){
    charge_He3 = 51.586;
    charge_C = 0.25;

    z0 = 1.26572;
  }
  
  double rho_N2 = 0.15; // in amagat
  double rho_N2_amagat = 0.15; // in amagat
  double rho_C = 2266;  // mg/cm3
  double rho_air = 1015 * 100; // Pa from JLab website

  double l_C = 0.01 * 2.54; // 0.01 inches converted to cm
  double l_cut = foil_cut * nfoil * 2.0 * 100; // foil cuts in cm, also used for he3

  double temp = 294; // in K
  double R = 8.3144; // From idea gas law

  rho_air /= R*temp*1e6;  // now in mol / cm^3
  double amagat_to_mol_cm3 = 44.615 / 1e6;
  double mol_to_mg_C = 28.02 * 1000;
  double mol_to_mg_air = 28.57 * 1000;

  rho_N2 *= amagat_to_mol_cm3 * mol_to_mg_C;
  rho_air *= mol_to_mg_air;

  //These are in units of mg / cm2
  double m2_N2 = rho_N2 * l_cut;
  double m2_C = ( rho_C * l_C * nfoil ) + ( rho_air * l_cut );

  TString jmgr_He3 = "../../config/" + cfg + "_He3.cfg";
  TString jmgr_optics = "../../config/" + cfg + "_optics.cfg";

  Utilities::KinConf kin_He3 = Utilities::LoadKinConfig(jmgr_He3,1);
  Utilities::KinConf kin_optics = Utilities::LoadKinConfig(jmgr_optics,1);
  analyzed_tree *T_He3 = Utilities::LoadAnalyzedRootFiles(kin_He3,1,0);
  analyzed_tree *T_optics = Utilities::LoadAnalyzedRootFiles(kin_optics,1,0);
  
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
  
  double N_C = 0;
  double N_C_acc = 0;
  double N_C_bg = 0;
  double N_C_bg_acc = 0;
  double N_1 = 0;
  double N_2 = 0;
  double N_3 = 0;
  double N_4 = 0;
  double N_QE = 0;
  double N_QE_acc = 0;
  double N_QE_bg = 0;
  double N_QE_bg_acc = 0;
  double N_QE1 = 0;
  double N_QE2 = 0;
  double N_QE3 = 0;
  double N_QE4 = 0;

  int nevent = 0;
  int maxevent = T_optics->fChain->GetEntries(); // Used to loop through the tree

  /////Set the histograms
  double hxmin = -4;
  double hxmax = 2.5;
  if(cfg == "GEN2"){
    hxmin = -6;
    hxmax = 3;
  }

  TH1F *hvz_optics = new TH1F("hvz_optics","Carbon Data;z_{tg} (m)",100,-0.3,0.3);
  TH1F *hvz_He3 = new TH1F("hvz_He3","He3 Data;z_{tg} (m)",100,-0.3,0.3);
  TH2F *hdxdy_optics = new TH2F("hdxdy_optics","Carbon HCal;#Deltay;#Deltax",150,-2,2,150,-6,6);
  TH2F *hdxdy_He3 = new TH2F("hdxdy_He3","He3 HCal;#Deltay;#Deltax",150,-2,2,150,-6,6);
  TH1F *hW2_optics = new TH1F("hW2_optics","Carbon W2;W^{2} (GeV)",100,-1,5);
  TH1F *hW2_He3 = new TH1F("hW2_He3","He3 W2;W^{2} (GeV)",100,-1,5);
  TH1F *hcoin_optics = new TH1F("hcoin_optics","Carbon Coincidence;Time (ns)",100,0,200);
  TH1F *hcoin_He3 = new TH1F("hcoin_He3","He3 Coincidence;Time (ns)",100,0,200);
  TH2F *hxysieve_optics = new TH2F("hxysieve_optics","Carbon Sieve;Sieve Y (m);Sieve X (m)",100,-0.15,0.15,100,-0.35,0.35);
  TH2F *hxysieve_He3 = new TH2F("hxysieve_He3","He3 Sieve;Sieve Y (m);Sieve X (m)",100,-0.15,0.15,100,-0.35,0.35);
  TH1F *hdx = new TH1F("hdx","",100,hxmin,hxmax);

  while(nevent < maxevent){
    T_optics->GetEntry(nevent++);  //Get data for one event

    double xsieve = T_optics->xtgt + T_optics->thtgt*z0;
    double ysieve = T_optics->ytgt + T_optics->phtgt*z0;

    ////// Define all the cuts we will use on the data  ////////////////
    bool good_W2 = T_optics->W2 > W2min && T_optics->W2 < W2max;
    bool good_dy_elas = T_optics->dy > dymin && T_optics->dy < dymax;
    bool good_dx_elas = T_optics->dx > dxmin && T_optics->dx < dxmax;
    bool good_coin_time = T_optics->coin_time > coin_min && T_optics->coin_time < coin_max;
    bool acc_coin_time = T_optics->coin_time > coin_bg_low && T_optics->coin_time < coin_bg_hi;
    bool dy_bg_cut = T_optics->dy < dy_bg_min || T_optics->dy > dy_bg_max;
    bool vz_cut = false;
    bool xysieve_cut = false;

    for(int ifoil = 0; ifoil < nfoil; ifoil++){
      bool vz_single_cut = T_optics->vz > foilz[ifoil] - foil_cut && T_optics->vz < foilz[ifoil] + foil_cut;
      vz_cut = vz_cut || vz_single_cut;
    }

    for(int ix=0; ix < nxsieve; ix++){
      for(int iy=0; iy < nysieve; iy++){
	bool xysieve_single_cut = pow(xsieve - xs_cent[ix],2) /(0.015*0.015) + pow(ysieve - ys_cent[iy],2) /(0.015*0.015)  < 1;
	xysieve_cut = xysieve_cut || xysieve_single_cut;
      }
    }

    //////////////////////////////////////////////////////////////////////
    
    hvz_optics->Fill(T_optics->vz);
    hW2_optics->Fill(T_optics->W2);
    hxysieve_optics->Fill(ysieve,xsieve);
    hcoin_optics->Fill(T_optics->coin_time);	

    N_1 += 1;
    if(good_W2 && vz_cut)
      N_2 += 1;
    if(good_coin_time)
      N_3 += 1;
    if(good_dy_elas && good_dx_elas)
      N_4 += 1;
    

    // All QE cuts
    if(good_W2 && vz_cut && xysieve_cut){

      if(good_dy_elas && good_dx_elas){
	if(acc_coin_time) N_C_acc++;
	if(good_coin_time) N_C++;
      }

      if(dy_bg_cut){
	if(acc_coin_time) N_C_bg_acc++;
	if(good_coin_time) N_C_bg++;
	
	hdx->Fill(T_optics->dx);
      }

      if(good_coin_time) hdxdy_optics->Fill(T_optics->dy,T_optics->dx);

    }

  }

  nevent = 0;
  maxevent = T_He3->fChain->GetEntries(); 

  while(nevent < maxevent){
    //while(nevent < 1000000){
    T_He3->GetEntry(nevent++);  //Get data for one event

    double xsieve = T_He3->xtgt + T_He3->thtgt*z0;
    double ysieve = T_He3->ytgt + T_He3->phtgt*z0;
    
    ////// Define all the cuts we will use on the data  ////////////////
    bool good_W2 = T_He3->W2 > W2min && T_He3->W2 < W2max;
    bool good_dy_elas = T_He3->dy > dymin && T_He3->dy < dymax;
    bool good_dx_elas = T_He3->dx > dxmin && T_He3->dx < dxmax;
    bool good_coin_time = T_He3->coin_time > coin_min && T_He3->coin_time < coin_max;
    bool acc_coin_time = T_He3->coin_time > coin_bg_low && T_He3->coin_time < coin_bg_hi;
    bool dy_bg_cut = T_He3->dy < dy_bg_min || T_He3->dy > dy_bg_max;
    bool vz_cut = false;
    bool xysieve_cut = false;

    for(int ifoil = 0; ifoil < nfoil; ifoil++){
      bool vz_single_cut = T_He3->vz > foilz[ifoil] - foil_cut && T_He3->vz < foilz[ifoil] + foil_cut;
      vz_cut = vz_cut || vz_single_cut;
    }

    for(int ix=0; ix < nxsieve; ix++){
      for(int iy=0; iy < nysieve; iy++){
	bool xysieve_single_cut = pow(xsieve - xs_cent[ix],2) /(0.015*0.015) + pow(ysieve - ys_cent[iy],2) /(0.015*0.015)  < 1;
	xysieve_cut = xysieve_cut || xysieve_single_cut;
      }
    }
    //////////////////////////////////////////////////////////////////////

    hvz_He3->Fill(T_He3->vz);
    hW2_He3->Fill(T_He3->W2);  
    hcoin_He3->Fill(T_He3->coin_time);  
    if(good_W2) hdxdy_He3->Fill(T_He3->dy,T_He3->dx);
    hxysieve_He3->Fill(ysieve,xsieve);

    N_QE1 += 1;
    if(good_W2 && vz_cut)
      N_QE2 += 1;
    if(good_coin_time)
      N_QE3 += 1;
    if(good_dy_elas && good_dx_elas)
      N_QE4 += 1;

    // QE Cuts
    if(good_W2 && vz_cut && xysieve_cut){

      if(good_dy_elas && good_dx_elas){
	if(acc_coin_time) N_QE_acc++;
	if(good_coin_time) N_QE++;
      }

      if(dy_bg_cut){
	if(acc_coin_time) N_QE_bg_acc++;
	if(good_coin_time) N_QE_bg++;
      }
      
      if(good_coin_time){
	hdxdy_He3->Fill(T_He3->dy,T_He3->dx);
      }
    }
    
  }

  

  double A = ( charge_He3 / charge_C ) * ( m2_N2 / m2_C ); //this is convienant
  double F = A * ( 1.0*(N_C - N_C_acc) / (N_QE - N_QE_acc) );
  double F_err = sqrt(A*A * ( N_C / (N_QE*N_QE) + N_C*N_C / (N_QE*N_QE*N_QE) ));

  double F_bg = A * ( 1.0*(N_C_bg - N_C_bg_acc) / (N_QE_bg - N_QE_bg_acc) );
  double F_bg_err = A*sqrt( (N_C + N_C_acc) / pow(N_QE - N_QE_acc,2) + pow(N_C - N_C_acc,2)*(N_QE + N_QE_acc) / pow(N_QE - N_QE_acc,4));

  DBInfo.NitrogenFraction = F;
  DBInfo.NitrogenFractionErr = F_err;
  
  DB_SetCorrections(DBInfo);
  
  cout<<"--------------- Nitrogen Info -----------"<<endl;
  cout<<"Air Length (cm) = "<<l_cut<<endl;
  cout<<"Carbon Length (cm) = "<<l_C<<endl;
  cout<<"Density N2 (amagat) = "<<rho_N2_amagat<<endl;
  cout<<"Density C (mg/cm^3) = "<<rho_C<<endl;
  cout<<"Density air (mg/cm^3) = "<<rho_air<<endl;
  cout<<"-----------------------------------------"<<endl;
  cout<<"Charge He3 (C) = "<<charge_He3<<endl;
  cout<<"Charge Carbon (C) = "<<charge_C<<endl;
  cout<<"N2 Mass Density (mg/cm2) = "<<m2_N2<<endl;
  cout<<"Carbon Mass Density (mg/cm2) = "<<m2_C<<endl;
  cout<<"-----------------------------------------"<<endl;
  cout<<"N Carbon = "<<N_C<<endl;
  cout<<"N Carbon Acc = "<<N_C_acc<<endl;
  cout<<"N He3 = "<<N_QE<<endl;
  cout<<"N He3 Acc = "<<N_QE_acc<<endl;
  cout<<"-----------------------------------------"<<endl;
  cout<<"Nitrogen Fraction = "<<F<<" &  Fraction Err = "<<F_err<<endl;


  ////// Make all the plots ///////////////////////////
  TLine *line1[nfoil];
  TLine *line2[nfoil];
  TBox *box[nfoil + 1];

  TCanvas *c1 = new TCanvas("c1","",1600,600);
  c1->Divide(2,1);
  
  c1->cd(1);
  hvz_optics->Draw();

  for(int ifoil = 0; ifoil < nfoil; ifoil++){
    double foil_min = foilz[ifoil] - foil_cut;
    double foil_max = foilz[ifoil] + foil_cut;
    
    DrawLines(c1->cd(1),foil_min,foil_max);
  }

  for(int ifoil = 0; ifoil < nfoil + 1; ifoil++){
    if(ifoil == 0) box[ifoil] = new TBox(c1->cd(1)->GetUxmin(), c1->cd(1)->GetUymin(), foilz[ifoil] - foil_cut, c1->cd(1)->GetUymax());
    else if(ifoil == nfoil) box[ifoil] = new TBox(foilz[ifoil - 1] + foil_cut, c1->cd(1)->GetUymin(), c1->cd(1)->GetUxmax(), c1->cd(1)->GetUymax());
    else box[ifoil] = new TBox(foilz[ifoil - 1] + foil_cut, c1->cd(1)->GetUymin(), foilz[ifoil] - foil_cut, c1->cd(1)->GetUymax());

    box[ifoil]->SetFillColorAlpha(kRed, 0.2);
    box[ifoil]->Draw("same");
  }

  TLegend *legend = new TLegend(0.61,0.81,0.88,0.88);
  legend->AddEntry(box[0],"Region Removed","f");
  legend->SetLineColor(0);
  legend->Draw("same");
  
  c1->cd(2);
  hvz_He3->Draw();

  TCanvas *c3 = new TCanvas("c3","",1600,600);
  c3->Divide(2,1);
  
  TEllipse *n_spot = new TEllipse(dy_n[0],dx_n[0],dy_n[1],dx_n[1]);
  n_spot->SetFillStyle(0);
  n_spot->SetLineWidth(3);
  n_spot->SetLineColor(kRed);

  c3->cd(1);
  hdxdy_He3->Draw("colz");
  n_spot->Draw("same");

  c3->cd(2);
  hdxdy_optics->Draw("colz");
  n_spot->Draw("same");

  TCanvas *c4 = new TCanvas("c4","",1600,600);
  c4->Divide(2,1);

  

  c4->cd(1);
  hW2_He3->Draw();
  DrawLines(c4->cd(1),W2min,W2max);

  c4->cd(2);
  hW2_optics->Draw();
  DrawLines(c4->cd(2),W2min,W2max);


  TCanvas *c5 = new TCanvas("c5","",1600,600);
  c5->Divide(2,1);
  
  c5->cd(1);
  hcoin_He3->Draw();
  DrawLines(c5->cd(1),coin_min,coin_max);

  c5->cd(2);
  hcoin_optics->Draw();
  DrawLines(c5->cd(2),coin_min,coin_max);


  TCanvas *c6 = new TCanvas("c6","",1600,600);
  c6->Divide(2,1);

  TEllipse *sieve_spot[nxsieve][nysieve];

  for(int ix=0; ix < nxsieve; ix++){
    for(int iy=0; iy < nysieve; iy++){
      sieve_spot[ix][iy] = new TEllipse(ys_cent[iy],xs_cent[ix],0.015,0.015);
      sieve_spot[ix][iy]->SetLineColor(kRed);
      sieve_spot[ix][iy]->SetFillStyle(0);
    }
  }
  
  c6->cd(1);
  hxysieve_optics->Draw("colz");

  for(int ix=0; ix < nxsieve; ix++){
    for(int iy=0; iy < nysieve; iy++){
      sieve_spot[ix][iy]->Draw("same");
    }
  }

  c6->cd(2);
  hxysieve_He3->Draw("colz");

  for(int ix=0; ix < nxsieve; ix++){
    for(int iy=0; iy < nysieve; iy++){
      sieve_spot[ix][iy]->Draw("same");
    }
  }

  TCanvas *c7 = new TCanvas("c7","",800,600);
  hdx->Draw();

  
  TString plot_dir = "../../plots/";
  TString plot_name = "Nitrogen_cont_" + cfg + ".pdf";

  TString outputfile = plot_dir + plot_name;

  c1->Print(outputfile + "(");
  c6->Print(outputfile);
  c3->Print(outputfile);
  c4->Print(outputfile);
  c5->Print(outputfile + ")");
  //c7->Print(outputfile + ")");
  

  //////////////////////////////////////////////////////////////
  
}
