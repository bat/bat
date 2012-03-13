#ifndef __BINOMIALMODEL__H
#define __BINOMIALMODEL__H

#include <BAT/BCModel.h>

// ---------------------------------------------------------
class BinomialModel : public BCModel
{
   public:

      // Constructors and destructor
      BinomialModel();
      BinomialModel(const char * name);
      ~BinomialModel();

      // Set total numer of events
      void SetNTotal(int n)
         { fNTotal = n; };

      // Set selected numer of events
      void SetNSelected(int n)
         { fNSelected = n; };

      // Methods to overload, see file BinomialModel.cxx
      void DefineParameters();
      double LogAPrioriProbability(const std::vector<double> &parameters);
      double LogLikelihood(const std::vector<double> &parameters);

   private:
      // the total number of events
      int fNTotal;

      // the selected (observed) number of events
      int fNSelected;
};
// ---------------------------------------------------------

#endif

