#ifndef __COMBINATIONMODEL__H
#define __COMBINATIONMODEL__H

#include <BAT/BCModel.h>

// ---------------------------------------------------------
class CombinationModel : public BCModel
{
 public:

  // Constructors and destructor
  CombinationModel();
  CombinationModel(const char * name);
  ~CombinationModel();

  // Methods to overload, see file CombinationModel.cxx
  void DefineParameters();
  double LogAPrioriProbability(std::vector <double> parameters);
  double LogLikelihood(std::vector <double> parameters);
};
// ---------------------------------------------------------

#endif

