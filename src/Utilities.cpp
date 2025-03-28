#include "../include/Utilities.h"

namespace Utilities {
  
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

  /* #################################################
     ##                HCAL Related                 ##  
     ################################################# */
  //_____________________________________
  TH2F *TH2FHCALface_rc(std::string name) {
    // returns TH2F for HCAL face (row,col)
    /* NOTE: HCAL block id (ibblk) starts from 1 and goes up to 288 but both
             HCAL row (rowblk) and column starts from 0 and goes up to 23 and 
             11, respectively. Extremely annoying! */
    TH2F *h = new TH2F(name.c_str(), ";HCAL columns;HCAL rows",
		       expconst::hcalcol, 0, expconst::hcalcol,
		       expconst::hcalrow, 0, expconst::hcalrow);
    return h;
  }
  //_____________________________________
  TH2F *TH2FHCALface_xy_data(std::string name) {
    // returns TH2F for HCAL face (x,y) [Data]
    double y_min = expconst::yHCAL_r_DB - expconst::hcalblk_w/2.;
    double y_max = expconst::yHCAL_l_DB + expconst::hcalblk_w/2.;
    double x_min = expconst::xHCAL_t_DB - expconst::hcalblk_h/2.;
    double x_max = expconst::xHCAL_b_DB + expconst::hcalblk_h/2.;
    TH2F *h = new TH2F(name.c_str(), ";y_{HCAL} (m);x_{HCAL} (m)",
		       expconst::hcalcol, y_min, y_max,
		       expconst::hcalrow, x_min, x_max);
    return h;
  }
  //_____________________________________
  TH2F *TH2FHCALface_xy_simu(std::string name) {
    // returns TH2F for HCAL face (x,y) [Simu]
    double y_min = expconst::yHCAL_r_DB_MC - expconst::hcalblk_w/2.;
    double y_max = expconst::yHCAL_l_DB_MC + expconst::hcalblk_w/2.;
    double x_min = expconst::xHCAL_t_DB_MC - expconst::hcalblk_h/2.;
    double x_max = expconst::xHCAL_b_DB_MC + expconst::hcalblk_h/2.;
    TH2F *h = new TH2F(name.c_str(), ";y_{HCAL} (m);x_{HCAL} (m)",
		       expconst::hcalcol, y_min, y_max,
		       expconst::hcalrow, x_min, x_max);
    return h;
  }
  //_____________________________________
  TH2F *TH2FdxdyHCAL(std::string name) {
    // returns TH2F for dxdyHCAL
    TH2F *h = new TH2F(name.c_str(), "; y_{HCAL} - y_{exp} (m); x_{HCAL} - x_{exp} (m)",
		       250, -2, 2, 250, -4, 6);
    return h;
  }
  //_____________________________________
  void DrawArea(vector<double> dimensions, int lcolor=2, int lwidth=4, int lstyle=9) {
    /* Draws four lines to represent a rectangular cut area */
    double top = dimensions[0];                 // -X axis
    double bottom = dimensions[1];              // +X axis
    double right = dimensions[2];               // -Y axis
    double left = dimensions[3];                // +Y axis
    TLine line;
    line.SetLineColor(lcolor); 
    line.SetLineWidth(lwidth); 
    line.SetLineStyle(lstyle);
    line.DrawLine(right, bottom, left, bottom); // bottom margin
    line.DrawLine(right, top, left, top);       // top margin
    line.DrawLine(right, top, right, bottom);   // right margin
    line.DrawLine(left, top, left, bottom);     // left margin
  }


  /* #################################################
     ##              Kinematic Histograms           ##  
     ################################################# */
  TH1F *TH1FhW(std::string name) {
    // returns W histogram
    TH1F *h = new TH1F(name.c_str(), "W Distribution (GeV)", 250,0,2);
    return h;
  }
  TH1F *TH1FhQ2(std::string name,       // Name of histogram
		TString conf) {             // SBS config
    // returns Q2 histogram
    int nbin=0; double hmin=-100, hmax=-100;
    if (conf=="GEN2" || conf=="GEN3" || conf=="GEN4" || conf=="GEN4b") { nbin=100; hmin=1.; hmax=4.; } 
    else if (conf=="GMN14") { nbin=100; hmin=5.; hmax=10.; }
    else cerr << "[Utilities::TH1FhQ2] Enter valid SBS config!!" << endl;
    TH1F *h = new TH1F(name.c_str(), "Q^{2} Distribution (GeV^{2})", 
		       nbin, hmin, hmax);
    return h;
  }

  TDatime SetTime(string time_str){

    TDatime time(time_str.c_str());
    time.Set(time.GetYear(),time.GetMonth(),time.GetDay(),time.GetHour(),time.GetMinute(),0);
      
    return time;
  }

  KinConf LoadKinConfig(TString config_file, bool is_data){
    
    KinConf kin_info;

    JSONManager *jmgr = new JSONManager(config_file.Data());

    //These variables are currently not in simulated files
    if(is_data){
      kin_info.IHWP_Flip = jmgr->GetValueFromKey<int>("IHWP_Flip");
      kin_info.dymin = jmgr->GetValueFromKey<double>("dymin");
      kin_info.dymax = jmgr->GetValueFromKey<double>("dymax");

      vector<double> coin_time_cut;
      jmgr->GetVectorFromKey<double>("coin_time", coin_time_cut);
      kin_info.coin_min = coin_time_cut[0];
      kin_info.coin_max = coin_time_cut[1];
      jmgr->GetVectorFromKey<int>("runnums",kin_info.runnums);
      kin_info.nruns = jmgr->GetValueFromKey<int>("Nruns_to_ana");
    }

    // seting up the desired SBS configuration
    kin_info.conf = jmgr->GetValueFromKey_str("GEN_config");
    kin_info.sbsmag = jmgr->GetValueFromKey<int>("SBS_magnet_percent");


    // Choosing the model of calculation
    // model 0 => uses reconstructed p as independent variable
    // model 1 => uses reconstructed angles as independent variable
    // model 2 => uses 4-vector calculation
    kin_info.model = jmgr->GetValueFromKey<int>("model");

    // choosing nucleon type 
    kin_info.Ntype = jmgr->GetValueFromKey_str("Ntype");

    // setting up global cuts
    std::string gcut = jmgr->GetValueFromKey_str("global_cut");
    kin_info.globalcut = gcut.c_str();
    //TTreeFormula *GlobalCut = new TTreeFormula("GlobalCut", globalcut, C);

    // HCAL cut definitions
    kin_info.sbs_kick = jmgr->GetValueFromKey<double>("sbs_kick");
    jmgr->GetVectorFromKey<double>("dx_p", kin_info.dx_p);
    jmgr->GetVectorFromKey<double>("dy_p", kin_info.dy_p);
    kin_info.Nsigma_dx_p = jmgr->GetValueFromKey<double>("Nsigma_cut_dx_p");
    kin_info.Nsigma_dy_p = jmgr->GetValueFromKey<double>("Nsigma_cut_dy_p");
    jmgr->GetVectorFromKey<double>("dx_n", kin_info.dx_n);
    jmgr->GetVectorFromKey<double>("dy_n", kin_info.dy_n);
    kin_info.Nsigma_dx_n = jmgr->GetValueFromKey<double>("Nsigma_cut_dx_n");
    kin_info.Nsigma_dy_n = jmgr->GetValueFromKey<double>("Nsigma_cut_dy_n");
    jmgr->GetVectorFromKey<double>("h_dxHCAL_lims", kin_info.hdx_lim);
    jmgr->GetVectorFromKey<double>("h_dyHCAL_lims", kin_info.hdy_lim);
    kin_info.hcal_voffset = jmgr->GetValueFromKey<double>("hcal_voffset");
    kin_info.hcal_hoffset = jmgr->GetValueFromKey<double>("hcal_hoffset");

    // Analysis cut limits
    kin_info.W2min = jmgr->GetValueFromKey<double>("W2min");
    kin_info.W2max = jmgr->GetValueFromKey<double>("W2max");
    
    
    // File information
    kin_info.rootfile_dir = jmgr->GetValueFromKey_str("rootfile_dir");
     

    return kin_info;

  }
  
  // is_data = 0/1 (data/simulation)
  TChain *LoadRawRootFiles(KinConf kin_info, bool is_data){
    
    TChain *T = new TChain("T");
    
    int nruns = kin_info.nruns;
    string rootfile_dir = kin_info.rootfile_dir;

    if(!is_data){
      TString rootfiles = rootfile_dir + "replayed*.root";
      T->Add(rootfiles);
    }
    else{
      vector<int> runnums = kin_info.runnums;
      if (nruns < 1 || nruns > runnums.size()) nruns = runnums.size();
      for (int i=0; i<nruns; i++) {
	string rfname = rootfile_dir + Form("/*%d*",runnums[i]);
	T->Add(rfname.c_str());
      }
    }
    
    return T;
  }


  analyzed_tree *LoadAnalyzedRootFiles(KinConf kin_info, bool is_data, bool is_reduced, bool is_pion = false){
    
    TChain *T = new TChain("Tout");
    
    TString file_dir = "/w/halla-scshelf2102/sbs/jeffas/GEN_analysis/scripts/outfiles/good_files/";
    if(is_pion) file_dir = "/lustre19/expphy/volatile/halla/sbs/jeffas/GEN_root/pion_test/";
    TString root_file;
    if(!is_data)
      root_file = "QE_sim_" + kin_info.conf + "_sbs" + kin_info.sbsmag + "p_nucleon_" + kin_info.Ntype + "_model" + kin_info.model + ".root";
    else
      root_file = "QE_data_" + kin_info.conf + "_sbs" + kin_info.sbsmag + "p_nucleon_" + kin_info.Ntype + "_model" + kin_info.model + ".root";

    T->Add(file_dir + root_file);

    analyzed_tree *tree_info = new analyzed_tree(T,is_data,is_reduced);
    
    return tree_info;
  }

}
