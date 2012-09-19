#ifndef __BCFBUBACKGROUND__H
#define __BCFBUBACKGROUND__H

// ---------------------------------------------------------

#include <string>

class TH1D;
class TF1;
class TH1;

// ---------------------------------------------------------
class BCFBUBackground
{
 public:

   /** \name Constructors and destructors */
   /** @{ */

	 /**
    * The default constructor. 
    * @param name The name of the background process. */
	BCFBUBackground(std::string name);

   /**
    * The default destructor. */
   ~BCFBUBackground();

   /** @} */

   /** \name Member functions (get) */
   /** @{ */

	 /**
		* @return The name of the background process. */
   std::string GetName()
      { return fName; };

   /**
    * @return The TH1D histogram. */
   TH1* GetHistogram()
      { return fHistogram; };

   /** @} */

   /** \name Member functions (set) */
   /** @{ */

   /**
    * Set the histogram.
    * @param hist The TH1 histogram. */
   void SetHistogram(TH1* hist);

   /** @} */

 private:

	 /** 
    * The TH1 histogram. */
   TH1* fHistogram;

   /**
    * The name of the background process. */
   std::string fName;

};
// ---------------------------------------------------------

#endif

