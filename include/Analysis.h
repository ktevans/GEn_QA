#ifndef Analysis_H
#define Analysis_H

// These are variables found in the analysis result
double u_n = -1.9130427;
double P_n = 0.96;    //Set for now, but check later
double P_p = -0.03;    //Set for now, but check later
double P_beam_avg = 0;
double P_beam_avg_err = 0;
double P_He3_avg = 0;
double P_He3_avg_err = 0;
double field_hor_ang = 0;
double field_vert_ang = 0;

// Values for expansion to calculate GE/GM
const int nexp = 6;
double T_avg[nexp] = {0};
double T_1_Q2_avg = 0;
double Q2_avg = 0;
double tau_avg = 0;
double epsilon_avg = 0;
double Px_avg = 0;
double Pz_avg = 0;
int count_avg_exp = 0;
double A_quad;
double B_quad;
double C_quad;
double GEGMp;
double GEGMp_err;
double GEGMn;
double GEGMn_stat_err;
double GEGMn_sys_err;
double GEn;
double GEn_stat_err;
double GEn_sys_err;
double GMn;
double GMn_err;

TVector3 TargetPolDirection;

// Basic function to get the integral from a histogram range
double GetYield(TH1F *h, double xmin, double xmax){
  int binmin = h->FindBin(xmin);
  int binmax = h->FindBin(xmax);
  return h->Integral(binmin, binmax);
}

// Sets the target polarization direction
void SetHe3PolAngle(){
  double phi_rad = field_vert_ang*TMath::Pi()/180.0;  
  double theta_rad = field_hor_ang*TMath::Pi()/180.0;
  
  TargetPolDirection.SetMagThetaPhi(1.0, theta_rad, phi_rad);
}

// Calcualtes the error assuming A = Np - Nm / Np + Nm
double CalcAsymErr(int Np, int Nm){
  double D = Np - Nm;
  double N = Np + Nm;
  return sqrt((N*N - D*D)/(N*N*N));
}

double CalcFractionErr(double Nx, double N){
  return sqrt(Nx/(N*N) + Nx*Nx/(N*N*N));
}


void UpdateAverageKinematics(analyzed_tree *T);
void GetGEGMFromTheory(bool is_neutron, double Q2, double &R, double &R_err);
void GetAFromQ2(bool is_neutron, double Q2, double &A, double &A_err);
double GetGEGMFromA(double A, double A_err);
void GetGEGMFromA_old(double A);

double GetAFromGEGM_theory(double GEGM,double etheta, double Q2);

// Class to handle the distributions and the fitting procedure
class distribution_fits {
 public:

  TString bg_shape_option;   // sets the type of background shape used
  bool bg_shape_set = false;

  // Histograms for all the fit distributions
  TH1F *hdx_data = NULL;
  TH1F *hdx_sim_p = NULL;
  TH1F *hdx_sim_n = NULL;
  TH1F *hdx_bg_fit;
  TH1F *hdx_bg_data = NULL;
  TH1F *hdx_total_fit;

  // Function sets the background shape from a string. 
  // Options: pol2, pol3, pol4, gaus, from data
  void SetBgShapeOption(TString bg_shape){
    bg_shape_option = bg_shape;
    bg_shape_set = true;
  };


  // Functions get the integrals from the fits of neutron/proton/bg distributinos
  int GetNYield(double xmin, double xmax){ return GetYield(hdx_sim_n,xmin,xmax);};
  int GetPYield(double xmin, double xmax){ return GetYield(hdx_sim_p,xmin,xmax);};
  int GetBgYield(double xmin, double xmax){ return GetYield(hdx_bg_fit,xmin,xmax);};

  // Set histogram shapes
  void SetDataShape(TH1F *h){hdx_data = h;};
  void SetPShape(TH1F *h){hdx_sim_p = h;};
  void SetNShape(TH1F *h){hdx_sim_n = h;};
  void SetBgShape(TH1F *h){hdx_bg_data = h;};
  
  // Get the histogram data
  TH1F *GetDataHist(){return hdx_data;};
  TH1F *GetPHist(){return hdx_sim_p;};
  TH1F *GetNHist(){return hdx_sim_n;};
  TH1F *GetBgHist(){return hdx_bg_fit;};
  TH1F *GetTotalHist(){return hdx_total_fit;};

  TString GetBgShape(){return bg_shape_option;};

  double fitbg_pol4(double *x, double *par);
  double fitsim(double *x, double *par);
  double fitbg_pol3(double *x, double *par);
  double fitbg_pol2(double *x, double *par);
  double fitbg_gaus(double *x, double *par);  
  
  void He3_fit_dists();

  // Constructor
  distribution_fits() { }
  
};

// Class to handle analyzing the asymmetry calculations
class analyzed_info {
 public:

  //This is used for counts for a single run
  struct run_info {
    int N_raw_p;
    int N_raw_m;
    double P_He3_sum;  // Used to avg for He3
    double P_He3_err_sum;  // Used to avg for He3
    double P_He3;
    double P_He3_err;
    double P_beam[2];
    int count_avg;
    bool run_init;
    
    //Starting values
    run_info() : N_raw_p(0), N_raw_m(0), P_He3_sum(0), P_He3_err_sum(0), P_He3(0), count_avg(0), run_init(false) {}
  };

  std::map<int, run_info> Asym_runs;

  void IterateRawCount(int key, int helicity){ // Add helicity raw counts
    if(helicity == 1) Asym_runs[key].N_raw_p++;
    else if(helicity == -1) Asym_runs[key].N_raw_m++;
  } 
  void UpdatePolHe3(int key, double Pol){
    Pol /= 100; //convert to decimal;
    double err = Pol*0.05;  // Hunter said to use 5% error for now
    Asym_runs[key].P_He3_sum += Pol;
    Asym_runs[key].P_He3_err_sum += err;
    Asym_runs[key].count_avg++;

    Asym_runs[key].P_He3 = Asym_runs[key].P_He3_sum / Asym_runs[key].count_avg;
    Asym_runs[key].P_He3_err = Asym_runs[key].P_He3_err_sum / Asym_runs[key].count_avg;
  }
  void SetAvgPol();
  void SetPolBeam(int key, TDatime *evtime, vector<vector<TDatime*>> PolTime, vector<vector<double>> PolVal);



  int N_raw_p = 0;
  int N_raw_m = 0;

  // accidental background +/- helicities
  int N_acc_p = 0;
  int N_acc_m = 0;
  int N_acc_cont = 0;

  // Pion +/- helicities
  int N_pion_p = 0;
  int N_pion_m = 0;
  int N_pion_cont = 0;

  // Number of protons
  int N_proton = 0;

  // Inelastic +/- helicities
  int N_in_p = 0;
  int N_in_m = 0;
  int N_in_cont = 0;
  
  // FSI +/- helicities
  int N_FSI_p = 0;
  int N_FSI_m = 0;

  // Asymmetry values for each contribution
  double A_raw;
  double A_acc = -1000;
  double A_p;
  double A_p_phys;
  double A_pion = -1000;
  double A_in = -1000;;
  double A_FSI;
  double A_phys;

  // Signal fractions for each contribution
  double f_acc = -1000;
  double f_N2 = -1000;
  double f_p;
  double f_pion = -1000;
  double f_in = -1000;
  double f_FSI;
  double f_n;

  // Asymmetry errors
  double A_raw_err;
  double A_acc_err = -1000;
  double A_in_err = -1000;
  double A_FSI_err;
  double A_p_err;
  double A_p_phys_err;
  double A_pion_err = -1000;
  double A_phys_stat_err;
  double A_phys_sys_err;
  double A_phys_tot_err;

  // Signal fraction errors
  double f_acc_err = -1000;
  double f_N2_err = -1000;
  double f_p_err;
  double f_pion_err = -1000;
  double f_in_err = -1000;
  double f_FSI_err;


  void SetNProton(int N){ 
    N_proton = N;
  }
  void IterateAccidentalCount(int helicity){ // Add helicity accidental background counts
    if(helicity == 1) N_acc_p++;
    else if(helicity == -1) N_acc_m++;
  } 
  void IterateAccidentalCont(){ // Add counts for accidentals inside QE cuts
    N_acc_cont++;
  } 
  
  void CalcAsymVals();
  
  double CalcAsymFractionErr(double A_f, double f_err){
    double P = P_He3_avg*P_beam_avg*P_n;
    return (P*A_phys - A_f) / (P*f_n)*f_err;
  }

  // Constructor
  analyzed_info() { }

};


#endif
