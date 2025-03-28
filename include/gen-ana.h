/* 
   gen-ana.h will include all the gen-ana libraries in the most
   efficient way. Just including this file in an analysis script
   will be enough to get access to all the gen-ana libraries.
   -----
   Sean Jeffas 11/4/2022 modified originaly from:
   P. Datta <pdbforce@jlab.org> Created 09-17-2022
*/

#include "Constants.h"                // namespace constant
#include "../src/JSONManager.cpp"
#include "../src/SetROOTVar.cpp"      // namespace setrootvar
#include "../src/TreeInfo.cpp"
#include "../src/Cut.cpp"             // namespace cut
#include "../src/Utilities.cpp"       // namespace util_pd
#include "../src/ExpConstants.cpp"    // namespace expconst & class SBSconfig
#include "../src/KinematicVar.cpp"    // namespace kine
#include "../src/DBparse.cpp"         // namespace DBparse
#include "../src/Analysis.cpp"
#include "../src/FF_param.cpp"



/* --- List of gmn-ana libraries --- */
// Cut.h          : namespace Cut
// Constants.h    : namespace constant
// KinematicVar.h : namespace kine
// ExpConstants.h : namespace expconst & class SBSconfig
// SetROOTVar.h   : namespace setrootvar
// Utilities.h    : namespace util_pd
