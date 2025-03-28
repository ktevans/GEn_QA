#include "../include/ExpConstants.h"


namespace expconst {

  double ebeam(TString config){
    if(config=="GMN1")
      return 1.916;
    else if(config=="GMN4")
      return 3.7278;
    else if(config=="GMN7")
      return 7.906;
    else if(config=="GMN11")
      return 9.91;
    else if(config=="GMN14" || config=="GMN8")
      return 5.965;
    else if(config=="GMN9")
      return 4.013;
    else if(config=="GEN2")
      return 4.291;
    else if(config=="GEN3")
      return 6.373;
    else if(config=="GEN4" || config=="GEN4b")
      return 8.448;
    else{
      std::cerr << "Enter a valid SBS configuration!" << std::endl;
      return -1;
    }
  }
  //--------------------------------------------
  double bbtheta(TString config){
    if(config=="GMN1")
      return 51.0;
    else if(config=="GMN4")
      return 36.0;
    else if(config=="GMN7")
      return 40.0;
    else if(config=="GMN11")
      return 42.0;
    else if(config=="GMN14")
      return 46.5;
    else if(config=="GMN8")
      return 26.5;
    else if(config=="GMN9")
      return 49.0;
    else if(config=="GEN2")
      return 29.5;
    else if(config=="GEN3")
      return 36.5;
    else if(config=="GEN4" || config=="GEN4b")
      return 35.0;
    else{
      std::cerr << "Enter a valid SBS configuration!" << std::endl;
      return -1;
    }
  }
  //--------------------------------------------
  double bbdist(TString config){
    if(config=="GMN1")
      return 1.8518;
    else if(config=="GMN4")
      return 1.7988;
    else if(config=="GMN7")
      return 1.84896;
    else if(config=="GMN11")
      return 1.55146;
    else if(config=="GMN14")
      return 1.84787;
    else if(config=="GMN8")
      return 1.97473;
    else if(config=="GMN9")
      return 1.550;
    else if(config=="GEN2" || config=="GEN3" || config=="GEN4" || config=="GEN4b")
      return 1.63;
    else{
      std::cerr << "Enter a valid SBS configuration!" << std::endl;
      return -1;
    }
  }
  //--------------------------------------------
  double sbstheta(TString config){
    if(config=="GMN1")
      return 33.5;
    else if(config=="GMN4")
      return 31.9;
    else if(config=="GMN7")
      return 16.1;
    else if(config=="GMN11")
      return 13.3;
    else if(config=="GMN14")
      return 17.3;
    else if(config=="GMN8")
      return 29.9;
    else if(config=="GMN9")
      return 22.5;
    else if(config=="GEN2")
      return 34.7;
    else if(config=="GEN3")
      return 22.1;
    else if(config=="GEN4" || config=="GEN4b")
      return 18.0;
    else{
      std::cerr << "Enter a valid SBS configuration!" << std::endl;
      return -1;
    }
  }
  //--------------------------------------------
  double sbsdist(TString config){
    if(config=="GMN1"||config=="GMN4"||config=="GMN7"||config=="GMN11"
       ||config=="GMN14"||config=="GMN8"||config=="GMN9")
      return 2.25;
    else if(config=="GEN2" || config == "GEN3" || config=="GEN4" || config=="GEN4b")
      return 2.8;
    else{
      std::cerr << "Enter a valid SBS configuration!" << std::endl;
      return -1;
    }
  }
  //--------------------------------------------
  double hcaldist(TString config){
    if(config=="GMN1")
      return 13.5;
    else if(config=="GMN4"||config=="GMN8"||config=="GMN9")
      return 11.0;
    else if(config=="GMN7"||config=="GMN14")
      return 14.0;
    else if(config=="GMN11")
      return 14.5;
    else if(config=="GEN2" || config=="GEN3" || config=="GEN4" || config=="GEN4b")
      return 17.0;
    else{
      std::cerr << "Enter a valid SBS configuration!" << std::endl;
      return -1;
    }
  }
  //--------------------------------------------
  double hcaltheta(TString config){
    if(config=="GMN1")
      return 33.5;
    else if(config=="GMN4")
      return 31.9;
    else if(config=="GMN7")
      return 16.1;
    else if(config=="GMN11")
      return 13.3;
    else if(config=="GMN14")
      return 17.3;
    else if(config=="GMN8")
      return 29.4;
    else if(config=="GMN9")
      return 22.0;
    else if(config=="GEN2")
      return 34.7;
    else if(config=="GEN3")
      return 21.6; 
    else if(config=="GEN4" || config=="GEN4b")
      return 18.0; 
    else{
      std::cerr << "Enter a valid SBS configuration!" << std::endl;
      return -1;
    }
  }

} //::expconst
