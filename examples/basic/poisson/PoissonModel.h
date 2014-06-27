#ifndef __POISSONMODEL__H
#define __POISSONMODEL__H

#include <BAT/BCModel.h>

// ---------------------------------------------------------
class PoissonModel : public BCModel
{
   public:

      // Constructors and destructor
      PoissonModel(const char * name);
      ~PoissonModel();

      // set number of observed events
      void SetNObs(unsigned nobs);

      // get number of observed events
      unsigned GetNObs() const
         { return fNObs; };

      // Method to overload, see file PoissonModel.cxx
      double LogLikelihood(const std::vector<double> & parameters);

   private:
      // number of observed events
      unsigned fNObs;
};
// ---------------------------------------------------------

#endif

