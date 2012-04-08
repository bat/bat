#ifndef __BCMTFCOMPARISONTOOL__H
#define __BCMTFCOMPARISONTOOL__H

#include <string>
#include <vector>

#include <TH1D.h>

// ---------------------------------------------------------
class BCMTFComparisonTool
{

   public:

      // Constructors and destructor
      BCMTFComparisonTool(const char * name);
      ~BCMTFComparisonTool();

      // setters


      // getters

      // return the name
      std::string GetName()
         { return fName; };

      // return number of contributions
      int GetNContributions()
         { return (int) fHistogramContainer.size(); };

      // misc
      // add a contribution
      void AddContribution(const char * name, TH1D hist);

      // add a contribution
      void AddContribution(const char * name, double centralvalue, double uncertainty);

      // print histograms to file
      void PrintHistograms(const char * filename);

      // draw overview
      void DrawOverview();

      // print overview to file
      void PrintOverview(const char * filename);

 private:

      // the name of the comparison
      std::string fName;

      // the name of the contributions
      std::vector<std::string> fNameContainer;

      // a container of histograms
      std::vector<TH1D *> fHistogramContainer;

      // a container of central values
      std::vector<double> fCentralValueContainer;

      // a container of uncertainties
      std::vector<double> fUncertaintyContainer;

};
// ---------------------------------------------------------

#endif

