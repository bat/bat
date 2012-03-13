#ifndef __POISSONMODEL__H
#define __POISSONMODEL__H

#include <BAT/BCModel.h>

// ---------------------------------------------------------
class PoissonModel : public BCModel
{
   public:

      // Constructors and destructor
      PoissonModel();
      PoissonModel(const char * name);
      ~PoissonModel();

      // set number of observed events
      int SetNObs(int nobs);

      // get number of observed events
      int GetNObs()
         { return fNObs; };

      // Methods to overload, see file PoissonModel.cxx
      void DefineParameters();
      double LogAPrioriProbability(const std::vector<double> &parameters);
      double LogLikelihood(const std::vector<double> &parameters);

   private:
      // number of observed events
      int fNObs;
};
// ---------------------------------------------------------

#endif

