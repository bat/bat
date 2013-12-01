#include <BAT/BCAux.h>
#include <BAT/BCH1D.h>
#include <BAT/BCH2D.h>

#include <TCanvas.h>

#include <iostream>

#include "BCPPDiagnostics.h"

int main()
{

   // set nicer style for drawing than the ROOT default
   BCAux::SetStyle();

   BCPPDiagnostics* p = new BCPPDiagnostics();

   p->OpenRootFile("output.root");

   p->PrintInfo();

   /*
   p->CalculateBatchQuantities();

   p->PrintBatchQuantities("batch.pdf");

   p->PrintLogProbability("logprob.pdf", "alllog");

   p->PrintTrajectory(0, "trajectory0.pdf", 1, 0, 10000);
   p->PrintTrajectory(1, "trajectory1.pdf", 1, 0, 10000);
   p->PrintTrajectory(0, 1, "trajectory01.pdf", 1, 0, 10000);
   */

   p->PrintAutocorrelation(0, "autocorr.pdf", 1, 20, -1);

   return 0;

}

