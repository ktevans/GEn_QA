#ifndef TREE_INFO_H
#define TREE_INFO_H


class analyzed_tree {
  public:
    TChain *fChain;
    
    int runnum;
    TDatime *datetime;
    //cuts
    bool WCut;
    bool pCut;
    bool nCut;
    bool fiduCut;
    bool coinCut;
    double ebeam;
    //kinematics
    double weight;
    double fnucl;
    double nu;
    double Q2;
    double W2;
    double dpel;
    double ephi;
    double etheta;
    double pcentral;
    //track
    double vz;
    double vx;
    double vy;
    double xtgt;
    double ytgt;
    double thtgt;
    double phtgt;
    double thetabend;
    double xfp;
    double yfp;
    double thfp;
    double phfp;
    double trP;
    double trPx;
    double trPy;
    double trPz;
    //BBCAL
    double ePS;
    double xPS;
    double eSH;
    double xSH;
    double ySH;
    //HCAL
    double eHCAL;
    double xHCAL;
    double yHCAL;
    double xHCAL_exp;
    double yHCAL_exp;
    double dx;
    double dy;
    //GRINCH
    double grinch_track;
    double grinch_clus_size;
    //Timing Information
    double coin_time;
    double hcal_time;
    double bbcal_time;
    int nhodo_clus;
    double hodo_time[1000];
    //BPM and Raster information
    double BPMAx;
    double BPMAy;
    double Rasterx;
    double Rastery;
    double Raster2x;
    double Raster2y;
    //Beam/Target information
    int helicity;
    int IHWP;
    double He3Pol;

    analyzed_tree(TChain *tree=0,bool is_data=1, bool is_reduced=1);
    //virtual Int_t    Cut(Long64_t entry);
    virtual Int_t    GetEntry(Long64_t entry);
    //virtual Long64_t LoadTree(Long64_t entry);
    virtual void     Init(TChain *tree,bool is_data, bool is_reduced);
    //virtual void     Loop();
    //virtual Bool_t   Notify();
    virtual void     Show(Long64_t entry = -1);
};


#endif
