#ifndef PD_UTIL_H
#define PD_UTIL_H

#include "TH1F.h"
#include "TH2F.h"
#include "TLine.h"
#include "TLatex.h"

#include "../include/ExpConstants.h"



namespace Utilities {

  struct KinConf {
    TString conf;
    int sbsmag;
    int model;
    int IHWP_Flip;
    TString Ntype;
    TCut globalcut;
    double sbs_kick;
    vector<double> dx_p;
    vector<double> dy_p;
    vector<double> dx_n;
    vector<double> dy_n;
    double Nsigma_dx_p;
    double Nsigma_dy_p;
    double Nsigma_dx_n;
    double Nsigma_dy_n;
    vector<double> hdx_lim;
    vector<double> hdy_lim;
    double hcal_voffset;
    double hcal_hoffset;
    double W2min;
    double W2max;
    double dymin;
    double dymax;
    double coin_min;
    double coin_max;
    string rootfile_dir;
    vector<int> runnums;
    int nruns;
  };


  /* #################################################
     ##                HCAL Histograms              ##  
     ################################################# */
  TH2F *TH2FHCALface_rc(std::string name);      // returns TH2F for HCAL face (row,col)
  TH2F *TH2FHCALface_xy_data(std::string name); // returns TH2F for HCAL face (x,y) [Data]
  TH2F *TH2FHCALface_xy_simu(std::string name); // returns TH2F for HCAL face (x,y) [Simu]
  TH2F *TH2FdxdyHCAL(std::string name);         // returns TH2F for dxdyHCAL

  // draws rectangular cut regions
  void DrawArea(vector<double> dimensions,      // a vector with extreme points
		int lcolor,  // Default = 2 
		int lwidth,  // Default = 4
		int lstyle); // Default = 9


  /* #################################################
     ##              Kinematic Histograms           ##  
     ################################################# */
  TH1F *TH1FhW(std::string name);   // returns W histogram
  TH1F *TH1FhQ2(std::string name,   // returns Q2 histogram
		int conf);   // SBS config

  void DrawLines(TVirtualPad *pad, double xmin, double xmax, Color_t color);
  TBox *DrawBox(TVirtualPad *pad, double xmin, double xmax, Color_t color);

  TDatime SetTime(string time_str);
  KinConf LoadKinConfig(TString config_file, bool is_data);
  TChain *LoadRawRootFiles(KinConf kin_info, bool is_data);
  analyzed_tree *LoadAnalyzedRootFiles(KinConf kin_info, bool is_data, bool is_reduced, bool is_pion);

}

#endif
