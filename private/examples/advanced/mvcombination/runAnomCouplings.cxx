#include <BAT/BCH1D.h>
#include <BAT/BCH2D.h>
#include <BAT/BCLog.h>
#include <BAT/BCAux.h>
#include <BAT/BCSummaryTool.h>

#include <TMatrixT.h>
#include <TH2D.h>
#include <TCanvas.h>

#include <iostream>
#include <fstream>

#include "AnomCouplings.h"

int main(int argc, char *argv[])
{

   // set nicer style for drawing than the ROOT default
   BCAux::SetStyle();

   // open log file
   BCLog::OpenLog("log.txt");
   BCLog::SetLogLevel(BCLog::detail);

   // create new AnomCouplings object
   AnomCouplings * m = new AnomCouplings();

   // set precision
   m->MCMCSetPrecision(BCIntegrate::kMedium);
   
   if (argc==2) {
     int isopen = m->ReadInput(argv[1]);
     if (!isopen) {
       std::cout << "Could not open file. Exit." << std::endl;
       return 1;
     }
   }
   else {
     std::cout << "No input file specified. Exit." << std::endl;
     return 1;
   }

	 m->MarginalizeAll();
	 m->FindMode( m->GetBestFitParameters() );
	 m->PrintAllMarginalized("marginalized_all.pdf");
	 m->PrintResults("results_all.txt");

	 m->GetSlice("gl", "gr", m->GetBestFitParameters(), 200)->Print("gl_gr.pdf", "BTcB3CS1gmode");

  // clean up
   delete m;
   
   // close log file
   BCLog::CloseLog();
   
   // no error
   return 0;
}

