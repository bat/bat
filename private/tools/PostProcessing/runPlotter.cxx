#include <BAT/BCAux.h>
#include <BAT/BCH1D.h>
#include <BAT/BCH2D.h>

#include <iostream>

#include "BCPPMarginalize.h"

int main()
{

   // set nicer style for drawing than the ROOT default
   BCAux::SetStyle();

   BCPPMarginalize* p = new BCPPMarginalize();

   p->OpenRootFile("output.root");

   p->PrintInfo();

   BCH1D* bhist1d = p->BuildMarginalized1D(0,100, -10, 10);
   BCH2D* bhist2d = p->BuildMarginalized2D(0,100, -10, 10, 1, 100, -10, 10);

   bhist1d->Print("test1d.pdf");
   bhist2d->Print("test2d.pdf");

   std::cout <<p->GetParameterValue(0, 0, 0) << std::endl;
   std::cout <<p->GetParameterValue(0, 0, 10) << std::endl;
   std::cout <<p->GetParameterValue(0, 0, 20) << std::endl;
   std::cout <<p->GetParameterValue(0, 0, 30) << std::endl;
   std::cout <<p->GetLogProbabilityValue(0, 120000-1) << std::endl;
   std::cout <<p->GetLogProbabilityValue(0, 100000-1) << std::endl;

   return 0;

}

