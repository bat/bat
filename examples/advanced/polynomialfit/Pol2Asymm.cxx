#include "Pol2Asymm.h"

#include <BAT/BCLog.h>
#include <BAT/BCParameter.h>

#include <TString.h>

// ---------------------------------------------------------
Pol2Asymm::Pol2Asymm(const char * name)
	: Pol1Asymm(name)
{
}

// ---------------------------------------------------------
Pol2Asymm::~Pol2Asymm()
{
}

// ---------------------------------------------------------
void Pol2Asymm::DefineParameters() {
   this -> AddParameter("p0",  0.  ,  2.2);    // index 0
   this -> AddParameter("p1", -0.1 ,  0.05);  // index 1
   this -> AddParameter("p2",  0.0 ,  0.001); // index 2

	 SetPriorConstantAll();

	 PrintSummary();
}
