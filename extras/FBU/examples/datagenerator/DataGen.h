#include <string>

class TH1D;
class TH2D;
class TF1;

// ---------------------------------------------------------
class DataGen
{
 public:

   /** \name Constructors and destructors */
   /** @{ */

	 /**
    * The default constructor. 
    * @param name The name of the project.
		* @param nbinstruth The number of truth bins.
		* @param truthmin The lower edge of the truth interval. 
		* @param truthmax The upper edge of the truth interval. 
		* @param nbinsreco The number of reco bins.
		* @param recomin The lower edge of the reco interval. 
		* @param recomax The upper edge of the reco interval. */
	DataGen(std::string name, int nbinstruth, double truthmin, double truthmax, int nbinsreco, double recomin, double recomax);

   /**
    * The default destructor. */
   ~DataGen();

   /** @} */

   /** \name Member functions (misc) */
   /** @{ */

	 /**
		* Fills all histograms
		* @param nevents The number of events for filling the migration matrix.
		* @param nsignal The number of signal events.
		* @param nbackground The number of background events.
		* @return An error code. */
	 int FillHistograms(int nevents, int nsignal, int nbackground=0);

	 /**
		* Write histograms to a ROOT file. 
		* @param filename The name of the file. 
		* @return An error code. */
	 int Write(std::string filename);

   /** @} */

   /** \name Member functions (get) */
   /** @{ */

	 /**
		* @return The name of the project. */
   std::string GetName()
		 { return fName; };

   /**
    * @return The TF1 function. */
   TF1* GetFuncTruthSignal()
	 { return fFuncTruthSignal; };

   /**
    * @return The TF1 function. */
   TF1* GetFuncBackground()
	 { return fFuncBackground; };

   /**
    * @return The TF1 function. */
   TF1* GetFuncResolution()
	 { return fFuncResolution; };

   /**
    * @return The TH1D histogram. */
   TH1D* GetHistData()
	 { return fHistData; };

   /**
    * @return The TH1D histogram. */
   TH1D* GetHistTruthSignal()
	 { return fHistTruthSignal; };

   /**
    * @return The TH1D histogram. */
   TH1D* GetHistBackground()
	 { return fHistBackground; };

   /**
    * @return The TH2D histogram. */
   TH1D* GetHistMigrationMatrix()
	 { return fHistMigrationMatrix; };
	 
   /**
    * @return The TH1D histogram. */
   TH1D* GetHistEfficiency()
	 { return fHistEfficiency; };

	 /**
		* @return The number of bins for the truth distribution. */
	 int GetNbinsTruth()
	 { return fNbinsTruth; };

	 /**
		* @return The lower edge of the truth interval. */
	 int GetTruthMin()
	 { return fTruthMin; };

	 /**
		* @return The upper edge of the truth interval. */
	 int GetTruthMax()
	 { return fTruthMax; };

	 /**
		* @return The number of bins for the reconstructed distribution. */
	 int GetNbinsReco()
	 { return fNbinsReco; };

	 /**
		* @return The lower edge of the reco interval. */
	 int GetRecoMin()
	 { return fRecoMin; };

	 /**
		* @return The upper edge of the reco interval. */
	 int GetRecoMax()
	 { return fRecoMax; };

   /** @} */

   /** \name Member functions (set) */
   /** @{ */

   /**
    * @param hist The TF1 function. */
   void SetFuncTruthSignal(TF1* func)
	 { fFuncTruthSignal = func; };

   /**
    * @param hist The TF1 function. */
   void SetFuncBackground(TF1* func)
	 { fFuncBackground = func; };

   /**
    * @param hist The TF1 function. */
   void SetFuncResolution(TF1* func)
	 { fFuncResolution = func; };

	 /**
		* @param The number of truth bins. */
	 void SetNbinsTruth(int nbins)
	 { fNbinsTruth = nbins; };

	 /**
		* @param truthmin The lower edge of the truth interval. 
		* @param truthmax The upper edge of the truth interval. */
	 void SetMinMaxTruth(double truthmin, double truthmax)
	 { fTruthMin = truthmin; 
		 fTruthMax = truthmax; };

	 /**
		* @param The number of reco bins. */
	 void SetNbinsReco(int nbins)
	 { fNbinsReco = nbins; };

	 /**
		* @param recomin The lower edge of the reco interval. 
		* @param recomax The upper edge of the reco interval. */
	 void SetMinMaxReco(double recomin, double recomax)
	 { fRecoMin = recomin; 
		 fRecoMax = recomax; };

   /** @} */

 private:

	 /**
		* Create all necessary histograms. 
		* @return An error code. */
	 int CreateHistograms();

   /**
    * The name of the project. */
   std::string fName;

	 /**
		* The pdf for generating truth signal. */
	 TF1* fFuncTruthSignal;

	 /**
		* The pdf for generating (already smeared) background. */
	 TF1* fFuncBackground;

	 /**
		* The resolution function used for smearing. */
	 TF1* fFuncResolution;

	 /**
		* The data histogram including smeared signal and background. */
	 TH1D* fHistData;

	 /**
		* The truth pdf for signal. */
	 TH1D* fHistTruthSignal;

	 /**
		* The truth pdf for signal. */
	 TH1D* fHistBackground;

	 /**
		* The migration matrix. */
	 TH2D* fHistMigrationMatrix;

	 /**
		* The efficiency. */
	 TH1D* fHistEfficiency; 

	 /**
		* The number of bins for the truth distribution. */
	 int fNbinsTruth;

	 /**
		* The lower edge of the truth interval. */
	 double fTruthMin;

	 /**
		* The upper edge of the truth interval. */
	 double fTruthMax;

	 /**
		* The number of bins for the reconstructed distriution. */
	 int fNbinsReco;

	 /**
		* The lower edge of the reco interval. */
	 double fRecoMin;

	 /**
		* The upper edge of the reco interval. */
	 double fRecoMax;

};
// ---------------------------------------------------------

