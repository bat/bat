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

  void SetFlagPhysicalConstraints(bool flag)
  { fFlagPhysicalConstraints = flag; }; 

  // BAT methods

  double LogLikelihood(const std::vector<double> &parameters);

  void MCMCIterationInterface();

 private:

  // flag for imposing physical constraints or not
  bool fFlagPhysicalConstraints; 

  // histogram containing posterior for FR
  TH1D* fHistFR;


};
// ---------------------------------------------------------

#endif

