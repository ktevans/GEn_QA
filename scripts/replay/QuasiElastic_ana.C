//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
//   Created by Provakar Datta
//   Modified by Sean Jeffas, sj9ry@virginia.edu
//     Last Modified August 7, 2024
//   Modified by Kate Evans, ktevans@jlab.org
//     Last Modified March 28, 2025
//   The purpose of this script is to take a configuraiton
//   file for some SBS experiment and to produced some
//   analyzed output with cuts on good elastic events.
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


DBparse::DBInfo DBInfo;

// Load database files
void getDB(TString cfg){

  cout<<"Attempting to load DB File"<<endl;
  cout<<"---------------------------------------------------------------"<<endl;

   vector<DBparse::DBrequest> request = {
     {"He3 Polarization","He3 target polarization", 1},
     {"Beam Polarization","Beam Polarization values",1},
     {"Helicity Quality","Helicity readback good? (0/1 = bad/good)",1},
     {"Moller Quality","Moller measurements known? (0/1 = no/yes)",1},
     {"Field Measurement","Magnetic field direction",1}
  };

  DBInfo.cfg = cfg;
  DBInfo.var_req = request;

  DB_load(DBInfo);

  cout<<"---------------------------------------------------------------"<<endl;

}


int QuasiElastic_ana(const std::string configfilename, std::string filebase="../outfiles/QE_data")
{

  //****************************************************************************
  // Read in config file and define basic variables

  string configdir = "../../config/";

  gErrorIgnoreLevel = kError; // Ignores all ROOT warnings

  // Define a clock to get macro processing time
  TStopwatch *sw = new TStopwatch(); sw->Start();

  // reading input config file ---------------------------------------
  Utilities::KinConf kin_info = Utilities::LoadKinConfig(configdir + configfilename,1);

  getDB(kin_info.conf);

  // parsing trees
  TChain *C = LoadRawRootFiles(kin_info, 1);

  // seting up the desired SBS configuration
  TString conf = kin_info.conf;
  int sbsmag = kin_info.sbsmag;
  SBSconfig sbsconf(conf, sbsmag);
  sbsconf.Print();

  //****************************************************************************

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

  double GEMpitch = 10.0*TMath::DegToRad();
  TVector3 GEMzaxis(-sin(GEMpitch),0,cos(GEMpitch));
  TVector3 GEMyaxis(0,1,0);
  TVector3 GEMxaxis = (GEMyaxis.Cross(GEMzaxis)).Unit();

  //****************************************************************************
  // setting up ROOT tree branch addresses

  // setting up global cuts
  TCut globalcut = kin_info.globalcut;
  TTreeFormula *GlobalCut = new TTreeFormula("GlobalCut", globalcut, C);

  std::cout << "Global cut: " + globalcut << std::endl;

  int maxNtr=1000;
  int maxhcal = 100;
  const int maxClus = 1000;

  C->SetBranchStatus("*",0); // Turn off all branches

  // Some CODA event information
  double evtime;
  setrootvar::setbranch(C,"g","evtime",&evtime);

  // bbcal sh clus var
  double eSH, xSH, ySH, atimeSH, idSH;
  std::vector<std::string> bbcalclvar = {"e","x","y","atimeblk","idblk"};
  std::vector<void*> bbcalclvar_mem = {&eSH,&xSH,&ySH,&atimeSH,&idSH};
  setrootvar::setbranch(C,"bb.sh",bbcalclvar,bbcalclvar_mem);

  // bbcal ps clus var
  double ePS, xPS, yPS, atimePS, idPS;
  std::vector<std::string> bbcalpsclvar = {"e","x","y","atimeblk","idblk"};
  std::vector<void*> bbcalpsclvar_mem = {&ePS,&xPS,&yPS,&atimePS,&idPS};
  setrootvar::setbranch(C,"bb.ps",bbcalpsclvar,bbcalpsclvar_mem);

  // hcal clus var
  double eHCAL[maxhcal], xHCAL[maxhcal], yHCAL[maxhcal], rblkHCAL[maxhcal], cblkHCAL[maxhcal], idblkHCAL[maxhcal],tdctimeHCAL[maxhcal],atimeHCAL[maxhcal];
  std::vector<std::string> hcalclvar = {"e","x","y","rowblk","colblk","idblk","tdctimeblk","atimeblk"};
  std::vector<void*> hcalclvar_mem = {&eHCAL,&xHCAL,&yHCAL,&rblkHCAL,&cblkHCAL,&idblkHCAL,&tdctimeHCAL,&atimeHCAL};
  setrootvar::setbranch(C, "sbs.hcal", hcalclvar, hcalclvar_mem);

  // grinch var
  double grinch_track, grinch_clus_size, grinch_hit_time[maxhcal], grinch_pmtnum[maxhcal], grinch_amp[maxhcal];
  std::vector<std::string> grinchvar = {"trackindex","size"};
  std::vector<void*> grinchvar_mem = {&grinch_track,&grinch_clus_size};
  setrootvar::setbranch(C, "bb.grinch_tdc.clus", grinchvar, grinchvar_mem);
  std::vector<std::string> grinchvar_hit = {"time","pmtnum","amp"};
  std::vector<void*> grinchvar_mem_hit = {&grinch_hit_time,&grinch_pmtnum,&grinch_amp};
  setrootvar::setbranch(C, "bb.grinch_tdc.hit", grinchvar_hit, grinchvar_mem_hit);
  int grinch_nhits;
  setrootvar::setbranch(C, "Ndata.bb.grinch_tdc.hit","time",&grinch_nhits);



  // hodoscope
  double hodo_time[maxClus];
  int nhodo_clus;
  setrootvar::setbranch(C,"bb.hodotdc.clus.bar.tdc","meantime",&hodo_time);
  setrootvar::setbranch(C,"Ndata.bb.hodotdc.clus.bar.tdc","meantime",&nhodo_clus);

  // track var
  double ntrack, p[maxNtr],px[maxNtr],py[maxNtr],pz[maxNtr],xTr[maxNtr],yTr[maxNtr],thTr[maxNtr],phTr[maxNtr],vx[maxNtr],vy[maxNtr],vz[maxNtr],xtgt[maxNtr],ytgt[maxNtr],thtgt[maxNtr],phtgt[maxNtr],xfp[maxNtr],yfp[maxNtr],thfp[maxNtr],phfp[maxNtr];
  std::vector<std::string> trvar = {"n","p","px","py","pz","vx","vy","vz","tg_x","tg_y","tg_th","tg_ph","r_x","r_y","r_th","r_ph"};
  std::vector<void*> trvar_mem = {&ntrack,&p,&px,&py,&pz,&vx,&vy,&vz,&xtgt,&ytgt,&thtgt,&phtgt,&xfp,&yfp,&thfp,&phfp};
  setrootvar::setbranch(C,"bb.tr",trvar,trvar_mem);

  // tdctrig variable (N/A for simulation)
  int tdcElemN;
  double tdcTrig[maxNtr], tdcElem[maxNtr];
  std::vector<std::string> tdcvar = {"tdcelemID","tdcelemID","tdc"};
  std::vector<void*> tdcvar_mem = {&tdcElem,&tdcElemN,&tdcTrig};
  setrootvar::setbranch(C,"bb.tdctrig",tdcvar,tdcvar_mem,1);

  //Beam helicity variables
  double helicity;
  setrootvar::setbranch(C,"scalhel","hel",&helicity);

  //IHWP State
  double IHWP;
  setrootvar::setbranch(C,"IGL1I00OD16_16","",&IHWP);

  //BPMA
  double BPMAx, BPMAy;
  std::vector<std::string> BPMAvar = {"x","y"};
  std::vector<void*> BPMAvar_mem = {&BPMAx,&BPMAy};
  setrootvar::setbranch(C,"Lrb.BPMA",BPMAvar,BPMAvar_mem);

  //Raster
  double rawcurx, rawcury;
  std::vector<std::string> Rastervar = {"x","y"};
  std::vector<void*> Rastervar_mem = {&rawcurx,&rawcury};
  setrootvar::setbranch(C,"Lrb.Raster.rawcur",Rastervar,Rastervar_mem);

  //Raster 2
  double rawcur2x, rawcur2y;
  std::vector<std::string> Raster2var = {"x","y"};
  std::vector<void*> Raster2var_mem = {&rawcur2x,&rawcur2y};
  setrootvar::setbranch(C,"Lrb.Raster2.rawcur",Raster2var,Raster2var_mem);

  // turning on the remaining branches we use for the globalcut
  C->SetBranchStatus("bb.gem.track.nhits", 1);
  C->SetBranchStatus("bb.etot_over_p", 1);
  C->SetBranchStatus("sbs.hcal.nclus", 1);

  //****************************************************************************
  // Define the output tree (histograms)

  // defining the outputfile
  TString outFile = Form("%s_" + sbsconf.GetSBSconf() + "_sbs%dp_nucleon_%s_model%d.root", filebase.c_str(),  sbsconf.GetSBSmag(), Ntype.Data(), model);
  TFile *fout = new TFile(outFile.Data(), "RECREATE");

  // defining histograms

  // Generic Histograms
  TH1F *h_W = Utilities::TH1FhW("h_W");
  TH1F *h_W_cut = Utilities::TH1FhW("h_W_cut");
  TH1F *h_W_acut = Utilities::TH1FhW("h_W_acut");
  TH1D *h_dpel = new TH1D("h_dpel",";p/p_{elastic}(#theta)-1;",100,-0.3,0.3);
  TH1F *h_Q2 = Utilities::TH1FhQ2("h_Q2", conf);
  vector<double> hdx_lim = kin_info.hdx_lim;
  vector<double> hdy_lim = kin_info.hdy_lim;
  TH1F *h_dxHCAL = new TH1F("h_dxHCAL","; x_{HCAL} - x_{exp} (m);",int(hdx_lim[0]),hdx_lim[1],hdx_lim[2]);
  TH1F *h_dyHCAL = new TH1F("h_dyHCAL","; y_{HCAL} - y_{exp} (m);",int(hdy_lim[0]),hdy_lim[1],hdy_lim[2]);

  // BBCal Histograms
  TH1F *h_pse = new TH1F("h_pse","; preshower energy [GeV];",150,0.,3.);
  TH1F *h_she = new TH1F("h_she","; shower energy [GeV];",200,0.,4.);
  TH1F *h_eovp = new TH1F("h_eovp","; BBCal best cluster (SH + PS) energy / track momentum;",100,0.5,1.5);

  // BB GEM Histograms

  // GRINCH Histograms
  TH2F *h_LE = new TH2F("h_LE", "GRINCH Leading Edge [ns]", 510,0,510,60,-30,30);

  // HCal Histograms
  TH2F *h2_rcHCAL = Utilities::TH2FHCALface_rc("h2_rcHCAL");
  TH2F *h2_dxdyHCAL = Utilities::TH2FdxdyHCAL("h2_dxdyHCAL");
  TH2F *h2_xyHCAL_p = Utilities::TH2FHCALface_xy_data("h2_xyHCAL_p");
  TH2F *h2_xyHCAL_n = Utilities::TH2FHCALface_xy_data("h2_xyHCAL_n");

  // Hodo Histograms
  TH1F *h_coin_time = new TH1F("h_coin_time", "Coincidence time (HCal - Hodo) (ns)", 200, 50, 250);
  TH1F *h_trigElem = new TH1F("h_trigElem", "TDC Trigger Element", 50, 0, 5);

  // SBS GEM Histograms

  // Target Histograms

  //****************************************************************************
  // Defining interesting ROOT tree branches
  TTree *Tout = new TTree("Tout", "");
  int T_runnum;         Tout->Branch("runnum", &T_runnum, "runnum/I");
  TDatime T_datetime;   Tout->Branch("datetime", "TDatime", &T_datetime);
  double T_ebeam;       Tout->Branch("ebeam", &T_ebeam, "ebeam/D");

  //cuts
  bool WCut;            Tout->Branch("WCut", &WCut, "WCut/B");
  bool pCut;            Tout->Branch("pCut", &pCut, "pCut/B");
  bool nCut;            Tout->Branch("nCut", &nCut, "nCut/B");
  bool fiduCut;         Tout->Branch("fiduCut", &fiduCut, "fiduCut/B");
  bool coinCut;         Tout->Branch("coinCut", &coinCut, "coinCut/B");

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
  double T_thetabend;   Tout->Branch("thetabend", &T_thetabend, "thetabend/D");
  double T_xfp;         Tout->Branch("xfp", &T_xfp, "xfp/D");
  double T_yfp;         Tout->Branch("yfp", &T_yfp, "yfp/D");
  double T_thfp;        Tout->Branch("thfp", &T_thfp, "thfp/D");
  double T_phfp;        Tout->Branch("phfp", &T_phfp, "phfp/D");
  double T_trP;         Tout->Branch("trP", &T_trP, "trP/D");
  double T_trPx;        Tout->Branch("trPx", &T_trPx, "trPx/D");
  double T_trPy;        Tout->Branch("trPy", &T_trPy, "trPy/D");
  double T_trPz;        Tout->Branch("trPz", &T_trPz, "trPz/D");

  //BBCAL
  double T_ePS;         Tout->Branch("ePS", &T_ePS, "ePS/D");
  double T_xPS;         Tout->Branch("xPS", &T_xPS, "xPS/D");
  double T_yPS;         Tout->Branch("yPS", &T_yPS, "yPS/D");
  double T_atimePS;     Tout->Branch("atimePS", &T_atimePS, "atimePS/D");
  double T_idPS;        Tout->Branch("idPS", &T_idPS, "idPS/D");
  double T_eSH;         Tout->Branch("eSH", &T_eSH, "eSH/D");
  double T_xSH;         Tout->Branch("xSH", &T_xSH, "xSH/D");
  double T_ySH;         Tout->Branch("ySH", &T_ySH, "ySH/D");
  double T_atimeSH;     Tout->Branch("atimeSH", &T_atimeSH, "atimeSH/D");
  double T_idSH;        Tout->Branch("idSH", &T_idSH, "idSH/D");

  //HCAL
  double T_eHCAL;       Tout->Branch("eHCAL", &T_eHCAL, "eHCAL/D");
  double T_xHCAL;       Tout->Branch("xHCAL", &T_xHCAL, "xHCAL/D");
  double T_yHCAL;       Tout->Branch("yHCAL", &T_yHCAL, "yHCAL/D");
  double T_xHCAL_exp;   Tout->Branch("xHCAL_exp", &T_xHCAL_exp, "xHCAL_exp/D");
  double T_yHCAL_exp;   Tout->Branch("yHCAL_exp", &T_yHCAL_exp, "yHCAL_exp/D");
  double T_dx;          Tout->Branch("dx", &T_dx, "dx/D");
  double T_dy;          Tout->Branch("dy", &T_dy, "dy/D");
  double T_theta_pq;    Tout->Branch("theta_pq", &T_theta_pq, "theta_pq/D");
  double T_idblkHCAL;   Tout->Branch("idblkHCAL", &T_idblkHCAL, "idblkHCAL/D");

  //GRINCH
  double T_grinch_track;              Tout->Branch("grinch_track", &T_grinch_track, "grinch_track/D");
  double T_grinch_clus_size;          Tout->Branch("grinch_clus_size", &T_grinch_clus_size, "grinch_clus_size/D");
  int T_grinch_nhits;                 Tout->Branch("grinch_nhits", &T_grinch_nhits, "grinch_nhits/I");
  double T_grinch_hit_time[maxClus];  Tout->Branch("grinch_hit_time", &T_grinch_hit_time, "grinch_hit_time[grinch_nhits]/D");
  double T_grinch_pmtnum[maxClus];    Tout->Branch("grinch_pmtnum", &T_grinch_pmtnum, "grinch_pmtnum[grinch_nhits]/D");
  double T_grinch_amp[maxClus];       Tout->Branch("grinch_amp", &T_grinch_amp, "grinch_amp[grinch_nhits]/D");

  //Timing Information
  double T_coin_time;            Tout->Branch("coin_time", &T_coin_time, "coin_time/D");
  double T_hcal_time;            Tout->Branch("hcal_time", &T_hcal_time, "hcal_time/D");
  //double T_bbcal_time;           Tout->Branch("bbcal_time", &T_bbcal_time, "bbcal_time/D");
  int T_nhodo_clus;              Tout->Branch("nhodo_clus", &T_nhodo_clus, "nhodo_clus/I");
  double T_hodo_time[maxClus];   Tout->Branch("hodo_time", &T_hodo_time, "hodo_time[nhodo_clus]/D");

  //BPM and Raster information
  double T_BPMAx;     Tout->Branch("BPMAx", &T_BPMAx, "BPMAx/D");
  double T_BPMAy;     Tout->Branch("BPMAy", &T_BPMAy, "BPMAy/D");
  double T_rawcurx;   Tout->Branch("Rasterx", &T_rawcurx, "Rasterx/D");
  double T_rawcury;   Tout->Branch("Rastery", &T_rawcury, "Rastery/D");
  double T_rawcur2x;  Tout->Branch("Raster2x", &T_rawcur2x, "Raster2x/D");
  double T_rawcur2y;  Tout->Branch("Raster2y", &T_rawcur2y, "Raster2y/D");

  //Beam/Target information
  int T_helicity;        Tout->Branch("helicity", &T_helicity, "helicity/I");
  int T_IHWP;            Tout->Branch("IHWP", &T_IHWP, "IHWP/I");
  double T_He3Pol;       Tout->Branch("He3Pol", &T_He3Pol, "He3Pol/D");

  //****************************************************************************
  // Define cuts and kinematic offsets

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
  double coin_min = kin_info.coin_min;
  double coin_max = kin_info.coin_max;
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

  //****************************************************************************
  // Loop through the tree

  std::cout << std::endl;
  long nevent = 0, nevents = C->GetEntries();
  int treenum = 0, currenttreenum = 0, currentrunnum = 0;
  int IHWP_run = -100;
  time_t run_time_unix;

  cout<<"Processing "<<nevents<<" events"<<endl;

  while (C->GetEntry(nevent++))
  {

    // print progress
    if( nevent % 1000 == 0 ) std::cout << nevent*100.0/nevents << "% \r";
    std::cout.flush();

    //**************************************************************************
    // apply global cuts efficiently (AJRP method)
    currenttreenum = C->GetTreeNumber();
    if (nevent == 1 || currenttreenum != treenum)
    {
      treenum = currenttreenum;
      GlobalCut->UpdateFormulaLeaves();

      //Get the run number
      string s = C->GetFile()->GetName();
      int start = s.find("_stream0");
      start -= 4;
      int end = start + 4;
      T_runnum = stoi(s.substr(start,end - start));

      //Get the time this run started
      auto* Run_Data = C->GetFile()->Get<THaRunBase>("Run_Data");
      TDatime run_time = Run_Data->GetDate();
      run_time.Set(run_time.GetYear(),run_time.GetMonth(),run_time.GetDay(),run_time.GetHour(),run_time.GetMinute(),0);
      run_time_unix = run_time.Convert();

      // We must loop over a small subset to get the correct IHWP state for the entire run
      //First we do some stuff to max sure we dont reach the file limit
      int file_nevents = C->GetTree()->GetEntries();
      int max_events = 8000;
      if(file_nevents < max_events) max_events = file_nevents - 100;

      int start_event = nevent;

      while (C->GetEntry(start_event++) && start_event < nevent + max_events)
      {
	       if(IHWP == 1 || IHWP == -1) IHWP_run = IHWP;
      }

      C->GetEntry(nevent - 1);

    } // END if nevent = 1 or currenttree = tree

    bool passedgCut = GlobalCut->EvalInstance(0) != 0;
    if (!passedgCut) continue;

    //**************************************************************************
    // Define some timing variables

    double bbcal_trig_time=0., hcal_trig_time=0.;

    for(int ihit=0; ihit<tdcElemN; ihit++)
    {
      if(tdcElem[ihit]==5) bbcal_trig_time=tdcTrig[ihit];
      if(tdcElem[ihit]==0) hcal_trig_time=tdcTrig[ihit];
      h_trigElem->Fill(tdcElem[ihit]);
    }

    //Calculate absolute time of this event
    double time_interval = 2.004;  //in ns, used to be 4, but i think the interval for GEn is 2.004ns
    int time_rel = evtime*time_interval*1e-9 / 60; //in min, rounded

    TDatime time_abs(run_time_unix + time_rel * 60); //Add the relative minutes

    //**************************************************************************
    // Start filling in the output tree variables

    auto it = DBInfo.He3Pol.find(time_abs); // Get Polarization from the table

    if(it == DBInfo.He3Pol.end())
    {
      T_He3Pol = -1;
    }
    else
    {
      T_He3Pol = it->second;
    }

    T_datetime = time_abs;

    //Timing Information
    T_hcal_time = atimeHCAL[0];
    T_idblkHCAL = idblkHCAL[0];
    //T_bbcal_time = atimeSH;

    T_nhodo_clus = nhodo_clus;
    for(int iclus = 0; iclus < nhodo_clus; iclus++)
    {
      T_hodo_time[iclus] = hodo_time[iclus];
    }

    double coin_time = atimeHCAL[0] - hodo_time[0];
    T_coin_time = coin_time;
    h_coin_time->Fill(coin_time);

    coinCut = coin_time > coin_min && coin_time < coin_max;

    //Grinch info
    T_grinch_track = grinch_track;
    T_grinch_clus_size = grinch_clus_size;
    T_grinch_nhits = grinch_nhits;
    for(int iclus = 0; iclus < grinch_nhits; iclus++)
    {
      T_grinch_hit_time[iclus] = grinch_hit_time[iclus];
      T_grinch_pmtnum[iclus] = grinch_pmtnum[iclus];
      T_grinch_amp[iclus] = grinch_amp[iclus];
    }

    //Beam helicity information
    T_helicity = helicity;
    T_IHWP = IHWP_run;

    //BPMs and Rasters
    T_BPMAx = BPMAx;
    T_BPMAy = BPMAy;
    T_rawcurx = rawcurx;
    T_rawcury = rawcury;
    T_rawcur2x = rawcur2x;
    T_rawcur2y = rawcur2y;

    //**************************************************************************
    // Calculate new tree variables

    // kinematic parameters
    double ebeam = sbsconf.GetEbeam();  // Expected beam energy (GeV) [Get it from EPICS, eventually]
    double ebeam_corr = ebeam; //- MeanEloss;
    double precon = p[0]; //+ MeanEloss_outgoing

    // constructing the 4 vectors
    // Reaction    : e + e' -> N + N'
    // Conservation: Pe + Peprime = PN + PNprime
    TVector3 vertex(0, 0, vz[0]);
    TLorentzVector Pe(0, 0, ebeam_corr, ebeam_corr);  // incoming e-
    TLorentzVector Peprime(px[0] * (precon/p[0]), py[0] * (precon/p[0]), pz[0] * (precon/p[0]), precon);  // scattered e-
    TLorentzVector PN;  // target nucleon [Ntype ??]
    kine::SetPN(Ntype, PN);

    double etheta = kine::etheta(Peprime);
    double ephi = kine::ephi(Peprime);
    double pcentral = kine::pcentral(ebeam_corr, etheta, Ntype.Data());

    double nu = 0.;  // energy of the virtual photon
    double pN_expect = 0.;  // expected recoil nucleon momentum
    double thetaN_expect = 0.;  // expected recoil nucleon theta
    double phiN_expect = ephi + constant::pi;

    // Different modes of calculation. Goal is to achieve the best resolution
    // model 0 = uses reconstructed p as independent variable
    // model 1 = uses reconstructed angles as independent variable
    // model 2 = uses 4-vector calculation
    TVector3 pNhat;  // 3-momentum of the recoil nucleon (Unit)
    double Q2recon, W2recon;

    if (model == 0)
    {
      nu = Pe.E() - Peprime.E();
      pN_expect = kine::pN_expect(nu, Ntype.Data());
      thetaN_expect = acos((Pe.E() - Peprime.Pz()) / pN_expect);
      pNhat = kine::qVect_unit(thetaN_expect, phiN_expect);
      Q2recon = kine::Q2(Pe.E(), Peprime.E(), etheta);
      W2recon = kine::W2(Pe.E(), Peprime.E(), Q2recon, Ntype.Data());
    }
    else if (model == 1)
    {
      nu = Pe.E() - pcentral;
      pN_expect = kine::pN_expect(nu, Ntype.Data());
      thetaN_expect = acos((Pe.E() - pcentral*cos(etheta)) / pN_expect);
      pNhat = kine::qVect_unit(thetaN_expect, phiN_expect);
      Q2recon = kine::Q2(Pe.E(), Peprime.E(), etheta);
      W2recon = kine::W2(Pe.E(), Peprime.E(), Q2recon, Ntype.Data());
    }
    else if (model == 2)
    {
      TLorentzVector q = Pe - Peprime; // 4-momentum of virtual photon
      nu = q.E();
      pNhat = q.Vect().Unit();
      Q2recon = -q.M2();
      W2recon = (PN + q).M2();
    }

    h_Q2->Fill(Q2recon);
    double Wrecon = sqrt(max(0., W2recon));
    double dpel = Peprime.E()/pcentral - 1.0;
    h_dpel->Fill(dpel);

    TVector3 enhat_tgt( thtgt[0], phtgt[0], 1.0 );
    enhat_tgt = enhat_tgt.Unit();

    TVector3 enhat_fp( thfp[0], phfp[0], 1.0 );
    enhat_fp = enhat_fp.Unit();

    TVector3 enhat_fp_rot = enhat_fp.X() * GEMxaxis + enhat_fp.Y() * GEMyaxis + enhat_fp.Z() * GEMzaxis;

    double thetabend = acos( enhat_fp_rot.Dot( enhat_tgt ) );

    T_ebeam = Pe.E();

    T_nu = nu;
    T_Q2 = Q2recon;
    T_W2 = W2recon;
    T_dpel = dpel;
    T_ephi = ephi;
    T_etheta = etheta;
    T_pcentral = pcentral;

    //**************************************************************************
    // Fill in more output tree variables

    T_vz = vz[0];
    T_vx = vx[0];
    T_vy = vy[0];
    T_xtgt = xtgt[0];
    T_ytgt = ytgt[0];
    T_thtgt = thtgt[0];
    T_phtgt = phtgt[0];
    T_thetabend = thetabend;
    T_xfp = xfp[0];
    T_yfp = yfp[0];
    T_thfp = thfp[0];
    T_phfp = phfp[0];
    T_trP = p[0];
    T_trPx = px[0];
    T_trPy = py[0];
    T_trPz = pz[0];

    T_ePS = ePS;
    T_xPS = xPS;
    T_yPS = yPS;
    T_atimePS = atimePS;
    T_idPS = idPS;
    T_eSH = eSH;
    T_xSH = xSH;
    T_ySH = ySH;
    T_atimeSH = atimeSH;
    T_idSH = idSH;

    T_eHCAL = eHCAL[0];
    T_xHCAL = xHCAL[0];
    T_yHCAL = yHCAL[0];

    // Expected position of the q vector at HCAL
    vector<double> xyHCAL_exp; // xyHCAL_exp[0] = xHCAL_exp & xyHCAL_exp[1] = yHCAL_exp
    double theta_pq;
    kine::GetxyHCALexpect(vertex, pNhat, HCAL_origin, HCAL_axes, xyHCAL_exp);
    double dx = xHCAL[0] - xyHCAL_exp[0];
    double dy = yHCAL[0] - xyHCAL_exp[1];

    T_xHCAL_exp = xyHCAL_exp[0];
    T_yHCAL_exp = xyHCAL_exp[1];
    T_dx = dx;
    T_dy = dy;

    //**************************************************************************
    // Define cuts

    // HCAL active area and safety margin cuts [Fiducial region]
    bool AR_cut = cut::inHCAL_activeA(xHCAL[0], yHCAL[0], hcal_active_area);
    bool FR_cut = cut::inHCAL_fiducial(xyHCAL_exp[0], xyHCAL_exp[1], sbs_kick, hcal_safety_margin);

    // HCAL cuts
    pCut = pow((dx-dx_p[0]) / (dx_p[1]*Nsigma_cut_dx_p), 2) + pow((dy-dy_p[0]) / (dy_p[1]*Nsigma_cut_dy_p), 2) <= 1.;
    nCut = pow((dx-dx_n[0]) / (dx_n[1]*Nsigma_cut_dx_n), 2) + pow((dy-dy_n[0]) / (dy_n[1]*Nsigma_cut_dy_n), 2) <= 1.;

    fiduCut = AR_cut && FR_cut;
    WCut = W2recon > W2min && W2recon < W2max;

    //**************************************************************************
    // Fill Histograms

    // W cut and coincidence cut
    if (WCut && coinCut)
    {
      // fiducial cut
      //if (fiduCut) { //No fiducial cut for GEN?
	     h_dxHCAL->Fill(dx);
	     h_dyHCAL->Fill(dy);
	     h2_rcHCAL->Fill(cblkHCAL[0], rblkHCAL[0]);
	     h2_dxdyHCAL->Fill(dy, dx);

	     if (pCut)
       {
         h2_xyHCAL_p->Fill(xyHCAL_exp[1], xyHCAL_exp[0] - sbs_kick);
       }

	     if (nCut)
       {
         h2_xyHCAL_n->Fill(xyHCAL_exp[1], xyHCAL_exp[0]);
       }
	     //}
    } // END W cut and coin cut

    h_LE->Fill(grinch_pmtnum[0],grinch_hit_time[0]);

    h_W->Fill(Wrecon);
    // fiducial cut but no W cut
    //if (fiduCut) { //No fiducial cut for GEN
    //h_W_cut->Fill(Wrecon);
    if (pCut||nCut)
    {
	     h_W_cut->Fill(Wrecon);

	     //if (WCut && coinCut)
       if (coinCut)
       {
         h_pse->Fill(ePS);
         h_she->Fill(eSH);
         h_eovp->Fill((eSH + ePS) / p[0]);
       }

    } // END spot cut
    else
    {
	     h_W_acut->Fill(Wrecon);
    }
      //}

    Tout->Fill();
  } // END while loop through tree

  std::cout << std::endl << std::endl;

  //****************************************************************************
  // Canvas 1: General Plots

  TCanvas *c1 = new TCanvas("c1", "General Plots", 1200, 1000);
  c1->Divide(2,2);

  c1->cd(1);
  h2_dxdyHCAL->Draw("colz");
  h2_dxdyHCAL->SetTitle("dx and dy [m] with W and coin cuts");
  TEllipse Ep_p;
  Ep_p.SetFillStyle(0);
  Ep_p.SetLineColor(2);
  Ep_p.SetLineWidth(2);
  Ep_p.DrawEllipse(dy_p[0], dx_p[0], Nsigma_cut_dy_p*dy_p[1], Nsigma_cut_dx_p*dx_p[1], 0,360,0);
  TEllipse Ep_n;
  Ep_n.SetFillStyle(0); Ep_n.SetLineColor(3); Ep_n.SetLineWidth(2);
  Ep_n.DrawEllipse(dy_n[0], dx_n[0], Nsigma_cut_dy_n*dy_n[1], Nsigma_cut_dx_n*dx_n[1], 0,360,0);

  c1->cd(2);
  h_W->Draw();
  h_W->SetTitle("W Distribution [GeV], Spot cuts (red) & Anti - Spot cuts (blue)");
  h_W->SetLineColor(1);
  h_W_cut->Draw("same");
  h_W_cut->SetLineColor(2);
  h_W_acut->Draw("same");

  c1->cd(3);
  h2_xyHCAL_p->Draw("colz");
  h2_xyHCAL_p->SetTitle("Protons on HCal - Active Area (Red) & Fiducial Region (Blue)");
  Utilities::DrawArea(hcal_active_area);
  Utilities::DrawArea(hcal_safety_margin,4);

  c1->cd(4);
  h2_xyHCAL_n->Draw("colz");
  h2_xyHCAL_n->SetTitle("Neutrons on HCal - Active Area (Red) & Fiducial Region (Blue)");
  Utilities::DrawArea(hcal_active_area);
  Utilities::DrawArea(hcal_safety_margin,4);

  c1->SaveAs("../../plots/GeneralPlots.pdf");

  //****************************************************************************
  // Canvas 2: BBCal Plots

  TCanvas *c2 = new TCanvas("c2", "BBCal Plots", 1200, 1000);
  c2->Divide(2,2);

  c2->cd(1);
  h_pse->Draw();
  h_pse->SetTitle("Preshower energy with global, W, coin, and spot cuts");

  c2->cd(2);
  h_she->Draw();
  h_she->SetTitle("Shower energy with global, W, coin, and spot cuts");

  c2->cd(3);
  h_eovp->Draw();
  h_eovp->SetTitle("E/p with global, W, coin, and spot cuts");

  c2->SaveAs("../../plots/BBCalPlots.pdf");

  //****************************************************************************
  // Canvas 3: Hodo Plots

  TCanvas *c3 = new TCanvas("c3", "Hodo Plots", 1200, 1000);
  c3->Divide(2,2);

  c3->cd(1);
  h_coin_time->Draw();
  c3->Update();
  TLine *coinMin = new TLine(coin_min, 0, coin_min, 100000000);
  coinMin->SetLineColor(kRed);
  coinMin->Draw();
  TLine *coinMax = new TLine(coin_max, 0, coin_max, 100000000);
  coinMax->SetLineColor(kRed);
  coinMax->Draw();
  gPad->Update();

  c3->cd(2);
  h_trigElem->Draw();
<<<<<<< HEAD

=======
  
>>>>>>> 3f21ccfe3ab4fa9b0b8503786df08cd3762574d9
  c3->SaveAs("../../plots/HodoPlots.pdf");

  //****************************************************************************
  // Canvas 4: GRINCH Plots

  TCanvas *c4 = new TCanvas("c4", "GRINCH Plots", 1200, 1000);
  c4->Divide(1,2);

  c4->cd(1);
  h_LE->Draw("colz");
  c4->SaveAs("../../plots/GRINCHPlots.pdf");

  //****************************************************************************
  // Send to output file

  cout << "------" << endl;
  cout << " Output file : " << outFile << endl;
  cout << "------" << endl << endl;

  sw->Stop();
  cout << "CPU time elapsed = " << sw->CpuTime() << " s. Real time = " << sw->RealTime() << " s. " << endl << endl;

  c1->Write();
  c2->Write();
  c3->Write();
  c4->Write();

  fout->Write();
  sw->Delete();

  return 0;
}
