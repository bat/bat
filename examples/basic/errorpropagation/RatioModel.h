#ifndef __RATIOMODEL__H
#define __RATIOMODEL__H

#include <BAT/BCModel.h>

// ---------------------------------------------------------
class RatioModel : public BCModel
{
   public:

      // Constructor and destructor
      RatioModel(const char * name);
      ~RatioModel();

      double LogLikelihood(const std::vector<double> &parameters);

   private:

	    double fRatio;						// holds ratio of the two model parameters.
};
// ---------------------------------------------------------

#endif

