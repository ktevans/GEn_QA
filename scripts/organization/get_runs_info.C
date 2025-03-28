#include "../../include/gen-ana.h"
#include "../../dflay/src/JSONManager.cxx"


void get_runs_info(const char *configfilename){

  JSONManager *jmgr = new JSONManager(configfilename);
  std::vector<int> runnums; jmgr->GetVectorFromKey<int>("runnums",runnums);

  gSystem->Exec("rm -f temp.txt");
  
  for(int irun=0; irun < runnums.size(); irun++){
    //for(int irun=0; irun < 10; irun++){
    gSystem->Exec(Form("ls /mss/halla/sbs/GEnII/raw/*%i* | sort -n | tail -1 >> temp.txt",runnums[irun]));
    cout<<irun*100.0/runnums.size()<<"%"<<endl;
  }

  ifstream tempfile("temp.txt");    
  string line;

  int seg_total = 0;
  
  while(getline(tempfile, line)){
    int stream = (int)line[line.length() - 3] - 48;
    int nseg = (int)line[line.length() - 1] - 48;

    seg_total += (nseg + 1) * (stream + 1);
  }

  cout<<endl;
  cout<<"total segemts: "<<seg_total<<endl;


  gSystem->Exec("rm -f temp.txt");

}
