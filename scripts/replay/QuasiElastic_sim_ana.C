//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
//   Created by Provakar Datta
//   Modified by Sean Jeffas, sj9ry@virginia.edu
//   Last Modified August 7, 2024
//
//
//   The purpose of this script analyze simulated data for
//   SBS experiments. It is almost identical to QuasiElastic_ana.C
//   Eventually the scripts shoudl probably be merged.
//
//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////

#include <vector>
#include <iostream>

#include "TCut.h"
#include "TH1F.h"
#include "TLatex.h"
#include "TChain.h"
#include "TVector3.h"
#include "TStopwatch.h"
#include "TTreeFormula.h"
#include "TLorentzVector.h"

#include "../../include/gen-ana.h"

int QuasiElastic_sim_ana(const std::string configfilename, std::string filebase="../outfiles/QE_sim")
{

  string configdir = "../../config/";
  gErrorIgnoreLevel = kError; // Ignores all ROOT warnings
  
  // Define a clock to get macro processing time
  TStopwatch *sw = new TStopwatch(); sw->Start();
  
  // reading input config file ---------------------------------------
  Utilities::KinConf kin_info = Utilities::LoadKinConfig(configdir + configfilename,0);
  
  // parsing trees
  TChain *C = LoadRawRootFiles(kin_info, 0);
  
  // seting up the desired SBS configuration
  TString conf = kin_info.conf;
  int sbsmag = kin_info.sbsmag;
  SBSconfig sbsconf(conf, sbsmag);
  sbsconf.Print();
 
  
  vector<int> Ntried;
  vector<double> luminosity, genvol;
  /*
  TString file = rootfile_dir + "simcout/"+conf+"_He3_"+SIMC_type+"_elastic_summary.csv";
  fstream file_csv; file_csv.open(file);
  //Check if file exists
  if (!file_csv.is_open()) {
    cout<<"!!!!WARNING, ERROR: File "<<file<<" is not found"<<endl;
    exit(0);
  }
  else {
    cout<<"DB file found: "<<file<<endl;
  }
  
  string line;
  int iline=0;
  while (getline(file_csv, line)) {
      // Create a stringstream to parse the line
    stringstream ss(line);
    string cell;
    iline++;   // Keep track of the number of lines
    if(iline == 1) continue; // Skip first line
    int irow = 0;
    // Split the line into cells using a comma as a delimiter
    while (getline(ss, cell, ',')) {
      irow++;
      if(irow == 3) Ntried.push_back(stoi(cell));
      if(irow == 4) genvol.push_back(stod(cell));
      if(irow == 5) luminosity.push_back(stod(cell));
    }
  }
  */
 
  double bbtheta = sbsconf.GetBBtheta_rad();

  TVector3 BB_zaxis( sin(bbtheta), 0.0, cos(bbtheta) );
  TVector3 BB_xaxis(0,-1,0);
  TVector3 BB_yaxis = (BB_zaxis.Cross(BB_xaxis)).Unit();
	  

  // Choosing the model of calculation
  // model 0 => uses reconstructed p as independent variable
  // model 1 => uses reconstructed angles as independent variable
  // model 2 => uses 4-vector calculation
  int model = kin_info.model;
  if (model == 0) std::cout << "Using model 0 [recon. p as indep. var.] for analysis.." << std::endl;
  else if (model == 1) std::cout << "Using model 1 [recon. angle as indep. var.] for analysis.." << std::endl;
  else if (model == 2) std::cout << "Using model 2 [4-vector calculation] for analysis.." << std::endl;
  else { std::cerr << "Enter a valid model number! **!**" << std::endl; throw; }
  
  
  // choosing nucleon type 
  TString Ntype = kin_info.Ntype;

  // setting up global cuts
  TCut globalcut = kin_info.globalcut;
  TTreeFormula *GlobalCut = new TTreeFormula("GlobalCut", globalcut, C);

  // setting up ROOT tree branch addresses ---------------------------------------
  int maxNtr=1000;
  C->SetBranchStatus("*",0);
  // beam energy - Probably we should take an average over 100 events
  // double HALLA_p;
  // setrootvar::setbranch(C, "HALLA_p", "", &HALLA_p);

  
  // bbcal clus var
  double eSH,atimeSH;
  std::vector<std::string> bbcalclvar = {"e","atimeblk"}; 
  std::vector<void*> bbcalclvar_mem = {&eSH,&atimeSH}; 
  setrootvar::setbranch(C,"bb.sh",bbcalclvar,bbcalclvar_mem);
  double ePS; setrootvar::setbranch(C,"bb.ps","e",&ePS);
  
  // hcal clus var
  double eHCAL, xHCAL, yHCAL, rblkHCAL, cblkHCAL, idblkHCAL,tdctimeHCAL;
  std::vector<std::string> hcalclvar = {"e","x","y","rowblk","colblk","idblk","tdctimeblk"};
  std::vector<void*> hcalclvar_mem = {&eHCAL,&xHCAL,&yHCAL,&rblkHCAL,&cblkHCAL,&idblkHCAL,&tdctimeHCAL};
  setrootvar::setbranch(C, "sbs.hcal", hcalclvar, hcalclvar_mem);
  
  // grinch var
  double grinch_track,grinch_clus_size;
  std::vector<std::string> grinchvar = {"trackindex","size"};
  std::vector<void*> grinchvar_mem = {&grinch_track,&grinch_clus_size};
  setrootvar::setbranch(C, "bb.grinch.clus", grinchvar, grinchvar_mem);

  // track var
  double ntrack, p[maxNtr],px[maxNtr],py[maxNtr],pz[maxNtr],xTr[maxNtr],yTr[maxNtr],thTr[maxNtr],phTr[maxNtr];
  double vx[maxNtr],vy[maxNtr],vz[maxNtr];
  double xtgt[maxNtr],ytgt[maxNtr],thtgt[maxNtr],phtgt[maxNtr];
  double xfp[maxNtr],yfp[maxNtr],thfp[maxNtr],phfp[maxNtr];
  std::vector<std::string> trvar = {"n","p","px","py","pz","vx","vy","vz","tg_x","tg_y","tg_th","tg_ph","r_x","r_y","r_th","r_ph"};
  std::vector<void*> trvar_mem = {&ntrack,&p,&px,&py,&pz,&vx,&vy,&vz,&xtgt,&ytgt,&thtgt,&phtgt,&xfp,&yfp,&thfp,&phfp};
  setrootvar::setbranch(C,"bb.tr",trvar,trvar_mem);

  //MC variables
  double mc_sigma, mc_omega, mc_fnucl, mc_simc_weight;
  std::vector<std::string> mc = {"mc_sigma","mc_omega","mc_fnucl","simc_Weight"};
  std::vector<void*> mc_mem = {&mc_sigma,&mc_omega,&mc_fnucl,&mc_simc_weight};
  setrootvar::setbranch(C,"MC",mc,mc_mem);

  double mc_vx[maxNtr], mc_vy[maxNtr], mc_px[maxNtr], mc_py[maxNtr], mc_pz[maxNtr];
  std::vector<std::string> mc_vertex = {"sdtrack_vx","sdtrack_vy","sdtrack_momx","sdtrack_momy","sdtrack_momz"};
  std::vector<void*> mc_vertex_mem = {&mc_vx,&mc_vy,&mc_px,&mc_py,&mc_pz};
  setrootvar::setbranch(C,"MC",mc_vertex,mc_vertex_mem);

  // turning on the remaining branches we use for the globalcut
  C->SetBranchStatus("bb.gem.track.nhits", 1);
  C->SetBranchStatus("bb.etot_over_p", 1);
  C->SetBranchStatus("sbs.hcal.nclus", 1);
  
  // defining the outputfile
  TString outFile = Form("%s_" + sbsconf.GetSBSconf() + "_sbs%dp_nucleon_%s_model%d.root", 
			 filebase.c_str(),  sbsconf.GetSBSmag(), Ntype.Data(), model);
  TFile *fout = new TFile(outFile.Data(), "RECREATE");

  // defining histograms
  TH1F *h_W = Utilities::TH1FhW("h_W");
  TH1F *h_W_cut = Utilities::TH1FhW("h_W_cut");
  TH1F *h_W_acut = Utilities::TH1FhW("h_W_acut");
  TH1D *h_dpel = new TH1D("h_dpel",";p/p_{elastic}(#theta)-1;",100,-0.3,0.3);
  
  TH1F *h_Q2 = Utilities::TH1FhQ2("h_Q2", conf);
  vector<double> hdx_lim = kin_info.hdx_lim;
  vector<double> hdy_lim = kin_info.hdy_lim;
  TH1F *h_dxHCAL = new TH1F("h_dxHCAL","; x_{HCAL} - x_{exp} (m);",int(hdx_lim[0]),hdx_lim[1],hdx_lim[2]);
  TH1F *h_dyHCAL = new TH1F("h_dyHCAL","; y_{HCAL} - y_{exp} (m);",int(hdy_lim[0]),hdy_lim[1],hdy_lim[2]);
  TH1F *h_coin_time = new TH1F("h_coin_time", "Coincidence time (ns)", 200, 380, 660);

  TH2F *h2_rcHCAL = Utilities::TH2FHCALface_rc("h2_rcHCAL");
  TH2F *h2_dxdyHCAL = Utilities::TH2FdxdyHCAL("h2_dxdyHCAL");
  TH2F *h2_xyHCAL_p = Utilities::TH2FHCALface_xy_data("h2_xyHCAL_p");
  TH2F *h2_xyHCAL_n = Utilities::TH2FHCALface_xy_data("h2_xyHCAL_n");

  // Defining interesting ROOT tree branches 
  TTree *Tout = new TTree("Tout", "");
  //cuts
  bool WCut;            Tout->Branch("WCut", &WCut, "WCut/B");
  bool pCut;            Tout->Branch("pCut", &pCut, "pCut/B");
  bool nCut;            Tout->Branch("nCut", &nCut, "nCut/B");
  double weight;        Tout->Branch("weight", &weight, "weight/D");
  double T_fnucl;         Tout->Branch("fnucl", &T_fnucl, "fnucl/D");
  bool fiduCut;         Tout->Branch("fiduCut", &fiduCut, "fiduCut/B");
  bool coinCut;         Tout->Branch("coinCut", &coinCut, "coinCut/B");
  //
  double T_ebeam;       Tout->Branch("ebeam", &T_ebeam, "ebeam/D");
  //kine
  double T_nu;          Tout->Branch("nu", &T_nu, "nu/D");
  double T_Q2;          Tout->Branch("Q2", &T_Q2, "Q2/D");
  double T_W2;          Tout->Branch("W2", &T_W2, "W2/D");
  double T_dpel;        Tout->Branch("dpel", &T_dpel, "dpel/D");
  double T_ephi;        Tout->Branch("ephi", &T_ephi, "ephi/D");
  double T_etheta;      Tout->Branch("etheta", &T_etheta, "etheta/D");
  double T_pcentral;    Tout->Branch("pcentral", &T_pcentral, "pcentral/D");
  //track
  double T_vz;          Tout->Branch("vz", &T_vz, "vz/D");
  double T_vx;          Tout->Branch("vx", &T_vx, "vx/D");
  double T_vy;          Tout->Branch("vy", &T_vy, "vy/D");
  double T_xtgt;        Tout->Branch("xtgt", &T_xtgt, "xtgt/D");
  double T_ytgt;        Tout->Branch("ytgt", &T_ytgt, "ytgt/D");
  double T_thtgt;       Tout->Branch("thtgt", &T_thtgt, "thtgt/D");
  double T_phtgt;       Tout->Branch("phtgt", &T_phtgt, "phtgt/D");
  double T_xfp;         Tout->Branch("xfp", &T_xfp, "xfp/D");
  double T_yfp;         Tout->Branch("yfp", &T_yfp, "yfp/D");
  double T_thfp;        Tout->Branch("thfp", &T_thfp, "thfp/D");
  double T_phfp;        Tout->Branch("phfp", &T_phfp, "phfp/D");  
  double T_trP;         Tout->Branch("trP", &T_trP, "trP/D");
  double T_trX;         Tout->Branch("trX", &T_trX, "trX/D");
  double T_trY;         Tout->Branch("trY", &T_trY, "trY/D");
  double T_trTh;        Tout->Branch("trTh", &T_trTh, "trTh/D");
  double T_trPh;        Tout->Branch("trPh", &T_trPh, "trPh/D");
  //BBCAL
  double T_ePS;         Tout->Branch("ePS", &T_ePS, "ePS/D"); 
  double T_eSH;         Tout->Branch("eSH", &T_eSH, "eSH/D"); 
  //HCAL
  double T_eHCAL;       Tout->Branch("eHCAL", &T_eHCAL, "eHCAL/D"); 
  double T_xHCAL;       Tout->Branch("xHCAL", &T_xHCAL, "xHCAL/D"); 
  double T_yHCAL;       Tout->Branch("yHCAL", &T_yHCAL, "yHCAL/D"); 
  double T_xHCAL_exp;   Tout->Branch("xHCAL_exp", &T_xHCAL_exp, "xHCAL_exp/D"); 
  double T_yHCAL_exp;   Tout->Branch("yHCAL_exp", &T_yHCAL_exp, "yHCAL_exp/D"); 
  double T_dx;          Tout->Branch("dx", &T_dx, "dx/D"); 
  double T_dy;          Tout->Branch("dy", &T_dy, "dy/D");
  //GRINCH
  double T_grinch_track;    Tout->Branch("grinch_track", &T_grinch_track, "grinch_track/D");
  double T_grinch_clus_size;    Tout->Branch("grinch_clus_size", &T_grinch_clus_size, "grinch_clus_size/D");

  // Do the energy loss calculation here ...........

  // HCAL cut definitions
  double sbs_kick = kin_info.sbs_kick;
  vector<double> dx_p = kin_info.dx_p;
  vector<double> dy_p = kin_info.dy_p;
  double Nsigma_cut_dx_p = kin_info.Nsigma_dx_p;
  double Nsigma_cut_dy_p = kin_info.Nsigma_dy_p;
  vector<double> dx_n = kin_info.dx_n;
  vector<double> dy_n = kin_info.dy_n;
  double Nsigma_cut_dx_n = kin_info.Nsigma_dx_n;
  double Nsigma_cut_dy_n = kin_info.Nsigma_dy_n;
  vector<double> coin_time_cut = kin_info.coin_time_cut;
  double Nsigma_coin_time = kin_info.Nsigma_coin_time;
  vector<double> hcal_active_area = cut::hcal_active_area_data(); // Exc. 1 blk from all 4 sides
  vector<double> hcal_safety_margin = cut::hcal_safety_margin(dx_p[1], dx_n[1], dy_p[1], hcal_active_area);

  // elastic cut limits
  double W2min = kin_info.W2min;
  double W2max = kin_info.W2max;

  // costruct axes of HCAL CoS in Hall CoS
  double hcal_voffset = kin_info.hcal_voffset;
  double hcal_hoffset = kin_info.hcal_hoffset;
  vector<TVector3> HCAL_axes; kine::SetHCALaxes(sbsconf.GetSBStheta_rad(), HCAL_axes);
  TVector3 HCAL_origin = sbsconf.GetHCALdist()*HCAL_axes[2] + hcal_voffset*HCAL_axes[0] + hcal_hoffset*HCAL_axes[1];

  // looping through the tree ---------------------------------------
  std::cout << std::endl;
  long nevent = 0, nevents = C->GetEntries(); 
  int treenum = 0, currenttreenum = 0;
  int run_inc = -1;
  
  cout<<"Processing "<<nevents<<" events"<<endl;
  
  while (C->GetEntry(nevent++)) {
    
    // print progress 
    if( nevent % 1000 == 0 ) std::cout << nevent*100.0/nevents << "% \r";
    std::cout.flush();
    
    // apply global cuts efficiently (AJRP method)
    currenttreenum = C->GetTreeNumber();
    if (nevent == 1 || currenttreenum != treenum) {
      treenum = currenttreenum;
      GlobalCut->UpdateFormulaLeaves();
      string s = C->GetFile()->GetName();
      int start = s.find("job_");
      start += 4;
      int end = s.find(".root");
      run_inc = stoi(s.substr(start,end - start));
    } 
    
    bool passedgCut = GlobalCut->EvalInstance(0) != 0;   
   
    if (!passedgCut) continue;
     
    //MC event weighting
    weight = mc_sigma*mc_omega;
    //weight = mc_simc_weight * luminosity[run_inc] * genvol[run_inc] * (1.0 / Ntried[run_inc]);
    
    //Grinch info
    T_grinch_track = grinch_track;
    T_grinch_clus_size = grinch_clus_size;

    // kinematic parameters
    double ebeam = sbsconf.GetEbeam();       // Expected beam energy (GeV) [Get it from EPICS, eventually]
    double ebeam_corr = ebeam; //- MeanEloss;
    double precon = p[0]; //+ MeanEloss_outgoing

    // constructing the 4 vectors
    // Reaction    : e + e' -> N + N'
    // Conservation: Pe + Peprime = PN + PNprime 
    TVector3 vertex(0, 0, vz[0]);
    TLorentzVector Pe(0,0,ebeam_corr,ebeam_corr);   // incoming e-
    TLorentzVector Peprime(px[0] * (precon/p[0]),   // scattered e-
			   py[0] * (precon/p[0]),
			   pz[0] * (precon/p[0]),
			   precon);                 
    TLorentzVector PN;                              // target nucleon [Ntype ??]
    kine::SetPN(Ntype.Data(), PN);

    double etheta = kine::etheta(Peprime);
    double ephi = kine::ephi(Peprime);
    double pcentral = kine::pcentral(ebeam_corr, etheta, Ntype.Data());

    double nu = 0.;                   // energy of the virtual photon
    double pN_expect = 0.;            // expected recoil nucleon momentum
    double thetaN_expect = 0.;        // expected recoil nucleon theta
    double phiN_expect = ephi + constant::pi; 
    // Different modes of calculation. Goal is to achieve the best resolution
    // model 0 = uses reconstructed p as independent variable
    // model 1 = uses reconstructed angles as independent variable 
    // model 2 = uses 4-vector calculation 
    TVector3 pNhat;                   // 3-momentum of the recoil nucleon (Unit)
    double Q2recon, W2recon;
    if (model == 0) {
      nu = Pe.E() - Peprime.E();
      pN_expect = kine::pN_expect(nu, Ntype.Data());
      thetaN_expect = acos((Pe.E() - Peprime.Pz()) / pN_expect);
      pNhat = kine::qVect_unit(thetaN_expect, phiN_expect);
      Q2recon = kine::Q2(Pe.E(), Peprime.E(), etheta);
      W2recon = kine::W2(Pe.E(), Peprime.E(), Q2recon, Ntype.Data());
    } else if (model == 1) {
      nu = Pe.E() - pcentral;
      pN_expect = kine::pN_expect(nu, Ntype.Data());
      thetaN_expect = acos((Pe.E() - pcentral*cos(etheta)) / pN_expect);
      pNhat = kine::qVect_unit(thetaN_expect, phiN_expect);
      Q2recon = kine::Q2(Pe.E(), Peprime.E(), etheta);
      W2recon = kine::W2(Pe.E(), Peprime.E(), Q2recon, Ntype.Data());
      /*
      cout<<nevent<<" "<<mc_px[0]<<" "<<px[0]<<endl;
      cout<<nevent<<" "<<mc_py[0]<<" "<<py[0]<<endl;
      cout<<nevent<<" "<<mc_pz[0]<<" "<<pz[0]<<endl;
      */
      
    } else if (model == 2) {
      TLorentzVector q = Pe - Peprime; // 4-momentum of virtual photon
      nu = q.E();
      pNhat = q.Vect().Unit();
      Q2recon = -q.M2();
      W2recon = (PN + q).M2();
    }
    h_Q2->Fill(Q2recon); 
    double Wrecon = sqrt(max(0., W2recon));
    double dpel = Peprime.E()/pcentral - 1.0; h_dpel->Fill(dpel);

    TVector3 pvect_true(mc_px[0],mc_py[0],mc_pz[0]);
    TVector3 pvect_BB( pvect_true.Dot(BB_xaxis), pvect_true.Dot(BB_yaxis), pvect_true.Dot(BB_zaxis) );
    
    
    T_ebeam = Pe.E();

    T_nu = nu;
    T_Q2 = Q2recon;
    T_W2 = W2recon;
    T_dpel = dpel;
    T_ephi = ephi;
    T_etheta = etheta;
    T_pcentral = pcentral;
    T_fnucl = mc_fnucl;

    T_vz = vz[0];
    T_vx = mc_vx[0];
    T_vy = mc_vy[0];
    T_xtgt = xtgt[0];
    T_ytgt = ytgt[0];
    T_thtgt = thtgt[0];
    T_phtgt = phtgt[0];
    T_xfp = xfp[0];
    T_yfp = yfp[0];
    T_thfp = thfp[0];
    T_phfp = phfp[0];
    T_trP = p[0];
    T_trX = xTr[0];
    T_trY = yTr[0];
    T_trTh = thTr[0];
    T_trPh = phTr[0];

    T_ePS = ePS;
    T_eSH = eSH;

    T_eHCAL = eHCAL;
    T_xHCAL = xHCAL;
    T_yHCAL = yHCAL;

    // Expected position of the q vector at HCAL
    vector<double> xyHCAL_exp; // xyHCAL_exp[0] = xHCAL_exp & xyHCAL_exp[1] = yHCAL_exp
    kine::GetxyHCALexpect(vertex, pNhat, HCAL_origin, HCAL_axes, xyHCAL_exp);
    double dx = xHCAL - xyHCAL_exp[0];  
    double dy = yHCAL - xyHCAL_exp[1]; 

    T_xHCAL_exp = xyHCAL_exp[0];
    T_yHCAL_exp = xyHCAL_exp[1];
    T_dx = dx;
    T_dy = dy;

    // HCAL active area and safety margin cuts [Fiducial region]
    bool AR_cut = cut::inHCAL_activeA(xHCAL, yHCAL, hcal_active_area);
    bool FR_cut = cut::inHCAL_fiducial(xyHCAL_exp[0], xyHCAL_exp[1], sbs_kick, hcal_safety_margin);

    // HCAL cuts
    pCut = pow((dx-dx_p[0]) / (dx_p[1]*Nsigma_cut_dx_p), 2) + pow((dy-dy_p[0]) / (dy_p[1]*Nsigma_cut_dy_p), 2) <= 1.;
    nCut = pow((dx-dx_n[0]) / (dx_n[1]*Nsigma_cut_dx_n), 2) + pow((dy-dy_n[0]) / (dy_n[1]*Nsigma_cut_dy_n), 2) <= 1.;

    fiduCut = AR_cut && FR_cut;
    WCut = W2recon > W2min && W2recon < W2max;

    // W cut and coincidence cut
    if (WCut && coinCut) {
      // fiducial cut
      //if (fiduCut) { //No fiducial cut for GEN?
      h_dxHCAL->Fill(dx,weight);
      h_dyHCAL->Fill(dy,weight);
	h2_rcHCAL->Fill(cblkHCAL, rblkHCAL);
	h2_dxdyHCAL->Fill(dy, dx);
	
	if (pCut) h2_xyHCAL_p->Fill(xyHCAL_exp[1], xyHCAL_exp[0] - sbs_kick);
	if (nCut) h2_xyHCAL_n->Fill(xyHCAL_exp[1], xyHCAL_exp[0]);
	//}
    }
    h_W->Fill(Wrecon);
    // fiducial cut but no W cut
    //if (fiduCut) { //No fiducial cut for GEN?
    //h_W_cut->Fill(Wrecon);
      if (pCut || nCut) { 
	h_W_cut->Fill(Wrecon);
      } else {
	h_W_acut->Fill(Wrecon);
      }
      //}
      
  
      Tout->Fill();
  } // event loop
  std::cout << std::endl << std::endl;
  
  TCanvas *c1 = new TCanvas("c1", "c1", 1200, 1000);
  c1->Divide(2,2);

  c1->cd(1); h2_dxdyHCAL->Draw("colz");
  TEllipse Ep_p;
  Ep_p.SetFillStyle(0); Ep_p.SetLineColor(2); Ep_p.SetLineWidth(2);
  Ep_p.DrawEllipse(dy_p[0], dx_p[0], Nsigma_cut_dy_p*dy_p[1], Nsigma_cut_dx_p*dx_p[1], 0,360,0);
  TEllipse Ep_n;
  Ep_n.SetFillStyle(0); Ep_n.SetLineColor(3); Ep_n.SetLineWidth(2);
  Ep_n.DrawEllipse(dy_n[0], dx_n[0], Nsigma_cut_dy_n*dy_n[1], Nsigma_cut_dx_n*dx_n[1], 0,360,0);
 
  c1->cd(2);
  h_W->Draw(); h_W->SetLineColor(1);
  h_W_cut->Draw("same"); h_W_cut->SetLineColor(2);
  h_W_acut->Draw("same");

  c1->cd(3);
  h2_xyHCAL_p->Draw("colz");
  Utilities::DrawArea(hcal_active_area);
  Utilities::DrawArea(hcal_safety_margin,4);

  c1->cd(4); 
  h2_xyHCAL_n->Draw("colz");
  Utilities::DrawArea(hcal_active_area);
  Utilities::DrawArea(hcal_safety_margin,4);

  // outFile.ReplaceAll(".root",".png");
  // c1->Print(outFile.Data(),"png");

  cout << "------" << endl;
  cout << " Output file : " << outFile << endl;
  cout << "------" << endl << endl;

  sw->Stop();
  cout << "CPU time elapsed = " << sw->CpuTime() << " s. Real time = " << sw->RealTime() << " s. " << endl << endl;

  c1->Write();
  
  fout->Write();
  sw->Delete();
    
  return 0;
}
