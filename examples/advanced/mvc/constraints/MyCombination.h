#ifndef __MYCOMBINATION__H
#define __MYCOMBINATION__H

#include <BAT/BCModel.h>
#include <BAT/BCMVCombination.h>

class TH1D;

// ---------------------------------------------------------
class MyCombination : public BCMVCombination
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

