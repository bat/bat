#ifndef __GAUSSMODEL__H
#define __GAUSSMODEL__H

#include <BAT/BCModel.h>

// ---------------------------------------------------------
class GaussModel : public BCModel
{
 public:

  // Constructors and destructor
  GaussModel();
  GaussModel(const char * name);
  ~GaussModel();

  // Methods to overload, see file GaussModel.cxx
  void DefineParameters();
  double LogAPrioriProbability(std::vector <double> parameters);
  double LogLikelihood(std::vector <double> parameters);
};
// ---------------------------------------------------------

#endif

