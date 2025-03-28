#ifndef DIST_H
#define DIST_H




class distribution_fits {
 public:

  TString bg_shape_option;
  bool bg_shape_set = false;

  TH1F *hdx_data = NULL;
  TH1F *hdx_sim_p = NULL;
  TH1F *hdx_sim_n = NULL;
  TH1F *hdx_bg_fit;
  TH1F *hdx_bg_data = NULL;
  TH1F *hdx_total_fit;

  void SetBgShapeOption(TString bg_shape){
    bg_shape_option = bg_shape;
    bg_shape_set = true;
  };
  
  
  void SetDataShape(TH1F *h){hdx_data = h;};
  void SetPShape(TH1F *h){hdx_sim_p = h;};
  void SetNShape(TH1F *h){hdx_sim_n = h;};
  void SetBgShape(TH1F *h){hdx_bg_data = h;};
  
  TH1F *GetDataHist(){return hdx_data;};
  TH1F *GetPHist(){return hdx_sim_p;};
  TH1F *GetNHist(){return hdx_sim_n;};
  TH1F *GetBgHist(){return hdx_bg_fit;};
  TH1F *GetTotalHist(){return hdx_total_fit;};

  TString GetBgShape(){return bg_shape_option;};
  double fitbg_pol4(double *x, double *par);
  double fitsim(double *x, double *par);
  double fitbg_pol3(double *x, double *par);
  double fitbg_gaus(double *x, double *par);  
  
  void He3_sim_fit();

  // Constructor
  distribution_fits() { }
  
};


#endif
