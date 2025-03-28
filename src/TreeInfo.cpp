#include "../include/TreeInfo.h"


analyzed_tree::analyzed_tree(TChain *tree,bool is_data, bool is_reduced) : fChain(0) 
{
  Init(tree, is_data, is_reduced);
}


// Inputs are like this:
// is_data = 0/1 (is data file / is simulation file)
// is_reduced = 0/1, reduced only reads the most important tree variables
void analyzed_tree::Init(TChain *tree,bool is_data, bool is_reduced)
{
  runnum = 0;
  datetime = new TDatime();
  WCut = false;
  pCut = false;
  nCut = false;
  fiduCut = false;
  coinCut = false;
  ebeam = 0;
  weight = 0;
  fnucl = 0;
  nu = 0;
  Q2 = 0;
  W2 = 0;
  dpel = 0;
  ephi = 0;
  etheta = 0;
  pcentral = 0;
  vz = 0;
  vx = 0;
  vy = 0;
  xtgt = 0;
  ytgt = 0;
  thtgt = 0;
  phtgt = 0;
  thetabend = 0;
  xfp = 0;
  yfp = 0;
  thfp = 0;
  phfp = 0;
  trP = 0;
  trPx = 0;
  trPy = 0;
  trPz = 0;
  ePS = 0;
  xPS = 0;
  eSH = 0;
  xSH = 0;
  ySH = 0;
  eHCAL = 0;
  xHCAL = 0;
  yHCAL = 0;
  xHCAL_exp = 0;
  yHCAL_exp = 0;
  dx = 0;
  dy = 0;
  grinch_track = 0;
  grinch_clus_size = 0;
  coin_time = 0;
  hcal_time = 0;
  bbcal_time = 0;
  nhodo_clus = 0;
  BPMAx = 0;
  BPMAy = 0;
  Rasterx = 0;
  Rastery = 0;
  Raster2x = 0;
  Raster2y = 0;
  helicity = 0;
  IHWP = 0;
  He3Pol = 0;

  if (!tree){
    cout<<"Error analyzed_tree::Init(): TTree is empty"<<endl;
    return;
  }

  fChain = tree;

  //fChain->SetBranchStatus("*",0);

  // These tree variables should always be read
  setrootvar::setbranch(fChain,"WCut","",&WCut);
  setrootvar::setbranch(fChain,"pCut","",&pCut);
  setrootvar::setbranch(fChain,"nCut","",&nCut);
  setrootvar::setbranch(fChain,"coinCut","",&coinCut);
  setrootvar::setbranch(fChain,"W2","",&W2);
  setrootvar::setbranch(fChain,"Q2","",&Q2);
  setrootvar::setbranch(fChain,"dx","",&dx);
  setrootvar::setbranch(fChain,"dy","",&dy);
  setrootvar::setbranch(fChain,"ePS","",&ePS);
  
  
  if(!is_data){ // These variables are only in simulation files
    setrootvar::setbranch(fChain,"weight","",&weight);
    setrootvar::setbranch(fChain,"fnucl","",&fnucl);
  }
  else { // These variables are only in the data files
    setrootvar::setbranch(fChain,"runnum","",&runnum);
    setrootvar::setbranch(fChain,"datetime","",&datetime);
    setrootvar::setbranch(fChain,"coin_time","",&coin_time);
    setrootvar::setbranch(fChain,"hcal_time","",&hcal_time);
    setrootvar::setbranch(fChain,"bbcal_time","",&bbcal_time);
    setrootvar::setbranch(fChain,"nhodo_clus","",&nhodo_clus);
    setrootvar::setbranch(fChain,"hodo_time","",&hodo_time);
    setrootvar::setbranch(fChain,"helicity","",&helicity);
    setrootvar::setbranch(fChain,"IHWP","",&IHWP);
    setrootvar::setbranch(fChain,"He3Pol","",&He3Pol);
  }
  if(!is_reduced) {// If not reduced than we read everything
    if(is_data){// Only data files has these
      setrootvar::setbranch(fChain,"thetabend","",&thetabend);
      setrootvar::setbranch(fChain,"trPx","",&trPx);
      setrootvar::setbranch(fChain,"trPy","",&trPy);
      setrootvar::setbranch(fChain,"trPz","",&trPz);
      setrootvar::setbranch(fChain,"xPS","",&xPS);
      setrootvar::setbranch(fChain,"xSH","",&xSH);
      setrootvar::setbranch(fChain,"ySH","",&ySH);
    }
    setrootvar::setbranch(fChain,"fiduCut","",&fiduCut);
    setrootvar::setbranch(fChain,"ebeam","",&ebeam);
    setrootvar::setbranch(fChain,"nu","",&nu);
    setrootvar::setbranch(fChain,"dpel","",&dpel);
    setrootvar::setbranch(fChain,"ephi","",&ephi);
    setrootvar::setbranch(fChain,"etheta","",&etheta);
    setrootvar::setbranch(fChain,"pcentral","",&pcentral);
    setrootvar::setbranch(fChain,"vz","",&vz);
    setrootvar::setbranch(fChain,"vx","",&vx);
    setrootvar::setbranch(fChain,"vy","",&vy);
    setrootvar::setbranch(fChain,"xtgt","",&xtgt);
    setrootvar::setbranch(fChain,"ytgt","",&ytgt);
    setrootvar::setbranch(fChain,"thtgt","",&thtgt);
    setrootvar::setbranch(fChain,"phtgt","",&phtgt);
    setrootvar::setbranch(fChain,"xfp","",&xfp);
    setrootvar::setbranch(fChain,"yfp","",&yfp);
    setrootvar::setbranch(fChain,"thfp","",&thfp);
    setrootvar::setbranch(fChain,"phfp","",&phfp);
    setrootvar::setbranch(fChain,"trP","",&trP);
    setrootvar::setbranch(fChain,"eSH","",&eSH);
    setrootvar::setbranch(fChain,"eHCAL","",&eHCAL);
    setrootvar::setbranch(fChain,"xHCAL","",&xHCAL);
    setrootvar::setbranch(fChain,"yHCAL","",&yHCAL);
    setrootvar::setbranch(fChain,"xHCAL_exp","",&xHCAL_exp);
    setrootvar::setbranch(fChain,"yHCAL_exp","",&yHCAL_exp);
    setrootvar::setbranch(fChain,"grinch_track","",&grinch_track);
    setrootvar::setbranch(fChain,"grinch_clus_size","",&grinch_clus_size);
  }
  
}


Int_t analyzed_tree::GetEntry(Long64_t entry)
{
  // Read contents of entry.
  if (!fChain) return 0;
  return fChain->GetEntry(entry);
}

void analyzed_tree::Show(Long64_t entry)
{
  // Print contents of entry.
  // If entry is not specified, print current entry
  if (!fChain) return;
  fChain->Show(entry);
}

