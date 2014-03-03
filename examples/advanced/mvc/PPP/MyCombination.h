#ifndef __MYCOMBINATION__H
#define __MYCOMBINATION__H

#include <BAT/BCModel.h>
#include <BAT/BCMVCombination.h>

class TH1D;
class TH2D;

// ---------------------------------------------------------
class MyCombination : public BCMVCombination
{
 public:

  // Constructor
  MyCombination();

  // Destructor
  ~MyCombination();

  // setters
  void SetHistRho(TH1D* hist)
  { fHistRho = hist; }; 

  void SetHistRhoAlpha(TH2D* hist)
  { fHistRhoAlpha = hist; }; 

  void SetHistRhoEta(TH2D* hist)
  { fHistRhoEta = hist; }; 

  void SetFlagPhysicalConstraints(bool flag)
  { fFlagPhysicalConstraints = flag; }; 

  // BAT methods

  double LogLikelihood(const std::vector<double> &parameters);

  void MCMCIterationInterface();

 private:

  // flag for imposing physical constraints or not
  bool fFlagPhysicalConstraints; 

  // histogram containing posterior for rho
  TH1D* fHistRho;
  TH2D* fHistRhoAlpha;
  TH2D* fHistRhoEta;

};
// ---------------------------------------------------------

#endif

