#ifndef __GAUSSMODEL__H
#define __GAUSSMODEL__H

#include <BAT/BCModel.h>

// ---------------------------------------------------------
class GaussModel : public BCModel
{
   public:

      // Constructor and destructor
      GaussModel(const char * name);
      ~GaussModel();

      // Method to overload, see file GaussModel.cxx
      double LogLikelihood(const std::vector<double> &parameters);
};
// ---------------------------------------------------------

#endif

