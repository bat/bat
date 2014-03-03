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
      double LogAPrioriProbability(const std::vector<double> &parameters);
      double LogLikelihood(const std::vector<double> &parameters);

      // overloaded trial function
      virtual double MCMCTrialFunctionSingle(unsigned int ichain, unsigned int ipar);
};
// ---------------------------------------------------------

#endif

