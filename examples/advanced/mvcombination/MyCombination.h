#ifndef __BAT__MyCOMBINATION__H
#define __BAT__MyCOMBINATION__H

#include <BAT/BCModel.h>

#include <TH1D.h>

#include "MVCombination.h"

// This is a MyCombination header file.
// Model source code is located in file MVCombination/MVCombination.cxx

// ---------------------------------------------------------
class MyCombination : public MVCombination
{
   public:

      // Constructor
      MyCombination();

			// Destructor
      ~MyCombination();

			// setters
			void SetHistFR(TH1D* hist)
			{ fHistFR = hist; }; 

			// BAT methods

      double LogAPrioriProbability(const std::vector<double> &parameters);

			void MCMCIterationInterface();

 private:

			TH1D* fHistFR;


};
// ---------------------------------------------------------

#endif

