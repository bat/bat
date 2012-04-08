#ifndef __BCMTFSYSTEMATICVARIATION__H
#define __BCMTFSYSTEMATICVARIATION__H

#include <string>
#include <vector>

class TH1D;

// ---------------------------------------------------------
class BCMTFSystematicVariation
{
 public:

   // Constructors and destructor
   BCMTFSystematicVariation(const char * channelname, const char * systematicname, int nprocesses);
   ~BCMTFSystematicVariation();

   // setters

   // set histogram
   void SetHistogramUp(int index, TH1D * hist)
      { fHistogramUpContainer[index] = hist; };

   // set histogram
   void SetHistogramDown(int index, TH1D * hist)
      { fHistogramDownContainer[index] = hist; };

   // set histogram
   void SetHistograms(int index, TH1D * hist_up, TH1D * hist_down)
      { fHistogramUpContainer[index] = hist_up;
        fHistogramDownContainer[index] = hist_down; };

   // getters

   // return histogram
   TH1D * GetHistogramUp(int index)
      { return fHistogramUpContainer.at(index); };

   // return histogram
   TH1D * GetHistogramDown(int index)
      { return fHistogramDownContainer.at(index); };

   // misc
   void AddHistogramUp(TH1D * hist)
      { fHistogramUpContainer.push_back(hist); };

   // misc
   void AddHistogramDown(TH1D * hist)
      { fHistogramDownContainer.push_back(hist); };

   // misc
   void AddHistograms(TH1D * hist_up, TH1D * hist_down)
      { fHistogramUpContainer.push_back(hist_up);
        fHistogramDownContainer.push_back(hist_down); };

 private:

   // a container of histograms
   std::vector<TH1D *> fHistogramUpContainer;

   // a container of histograms
   std::vector<TH1D *> fHistogramDownContainer;

   // channel name
   std::string fChannelName;

   // systematic name
   std::string fSystematicName;

};
// ---------------------------------------------------------

#endif

