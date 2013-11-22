#include <BAT/BCAux.h>
#include <BAT/BCH1D.h>
#include <BAT/BCH2D.h>

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

   p->CloseRootFile();

   return 0;

}

