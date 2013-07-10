#ifndef __RATIOMODEL__H
#define __RATIOMODEL__H

#include <BAT/BCModel.h>
#include <BAT/BCH1D.h>

// ---------------------------------------------------------
class RatioModel : public BCModel
{
   public:

      // Constructors and destructor
      RatioModel();
      RatioModel(const char * name);
      ~RatioModel();

      // Methods to overload, see file RatioModel.cxx
      void DefineParameters();
      void DefineHistogram();
      void PrintHistogram();
      double LogAPrioriProbability(const std::vector<double> &parameters);
      double LogLikelihood(const std::vector<double> &parameters);
      void MCMCUserIterationInterface();

   private:
      BCH1D * fHistRatio;
};
// ---------------------------------------------------------

#endif

