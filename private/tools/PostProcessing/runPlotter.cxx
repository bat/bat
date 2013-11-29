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

   p->CalculateBatchQuantities();

   p->PrintBatchQuantities("batch.pdf");

   p->PrintLogProbability("logprob.pdf", "all");

   return 1;

   return 0;

}

