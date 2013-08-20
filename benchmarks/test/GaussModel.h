
#ifndef __BC_TEST__GAUSSMODEL__H
#define __BC_TEST__GAUSSMODEL__H

#include <BAT/BCModel.h>

// ---------------------------------------------------------
class GaussModel : public BCModel
{
   public:

      // Constructors and destructor
      GaussModel(const char * name, const unsigned & nParameters, long loopIterations = 0l);
      virtual ~GaussModel();

      // Methods to overload, see file GaussModel.cxx
      double LogLikelihood(const std::vector<double> & parameters);
private:
      /**
       * Used in likelihood to prolong artificially.
       */
      unsigned long fLoopIterations;
};
// ---------------------------------------------------------

#endif
