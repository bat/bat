#ifndef __BCMTFCHANNEL__H
#define __BCMTFCHANNEL__H

#include <string>
#include <vector>

class BCMTFTemplate;
class BCMTFSystematicVariation;

// ---------------------------------------------------------
class BCMTFChannel
{

   public:

      // Constructors and destructor
      BCMTFChannel(const char * name);
      ~BCMTFChannel();

      // setters

      // set name
      void SetName(const char * name)
         { fName = name; };

      // set data
      void SetData(BCMTFTemplate * bctemplate)
         { fData = bctemplate; };

      // set a flag for using this channel or not
      void SetFlagChannelActive(bool flag)
         { fFlagChannelActive = flag; };

      // getters
      std::string GetName()
         { return fName; };

      // return a template
      BCMTFTemplate * GetTemplate(int index)
         { return fTemplateContainer.at(index); };

      // return a systematicvariation
      BCMTFSystematicVariation * GetSystematicVariation(int index)
         { return fSystematicVariationContainer.at(index); };

      // return the data
      BCMTFTemplate * GetData()
         { return fData; };

      // return flag
      bool GetFlagChannelActive()
         { return fFlagChannelActive; };

      // misc

      // add a template
      void AddTemplate(BCMTFTemplate * bctemplate)
         { fTemplateContainer.push_back(bctemplate); };

      // add a systematic variation
      void AddSystematicVariation(BCMTFSystematicVariation * variation)
         { fSystematicVariationContainer.push_back(variation); };

      // print templates
      void PrintTemplates(const char * filename);

      // print template with systematics
      void PrintTemplate(int index, const char * filename);

 private:

      // name of the channel
      std::string fName;

      // the data set
      BCMTFTemplate * fData;

      // a container of templates
      std::vector<BCMTFTemplate *> fTemplateContainer;

      // a container of systematics
      std::vector<BCMTFSystematicVariation *> fSystematicVariationContainer;

      // flag: channel is used (true) or not (false) in fit
      bool fFlagChannelActive;

};
// ---------------------------------------------------------

#endif

