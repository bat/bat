// ***************************************************************
// This file was created using the |:PROGRAM:| script
// for project |:Project:|.
// |:PROGRAM:| is part of Bayesian Analysis Toolkit (BAT).
// BAT can be downloaded from http://mpp.mpg.de/bat
// ***************************************************************

#include <BAT/BCLog.h>
#include <BAT/BCAux.h>

int main()
{
   // set nicer style for drawing than the ROOT default
   BCAux::SetStyle();

   // open log file
   BCLog::OpenLog("log.txt", BCLog::detail, BCLog::detail);

   // perform your analysis here

   // close log file
   BCLog::CloseLog();

   BCLog::OutSummary("Exiting");

   return 0;
}
