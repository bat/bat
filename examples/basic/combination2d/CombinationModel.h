#ifndef __COMBINATIONMODEL__H
#define __COMBINATIONMODEL__H

#include <BAT/BCModel.h>

// ---------------------------------------------------------
class CombinationModel : public BCModel
{
   public:

      // Constructor and destructor
      CombinationModel(const char * name);
      ~CombinationModel();

      // Method to overload, see file CombinationModel.cxx
      double LogLikelihood(const std::vector<double> &parameters);
};
// ---------------------------------------------------------

#endif

