#ifndef __PRIORMODEL__H
#define __PRIORMODEL__H

#include <BAT/BCModel.h>

// ---------------------------------------------------------
class PriorModel : public BCModel
{
	public:

		// Constructors and destructor
		PriorModel();
		PriorModel(const char * name);
		~PriorModel();

		void SetTestModel(BCModel* model)
		{ fTestModel = model; }; 

		int PerformAnalysis(); 

		double LogAPrioriProbability(std::vector <double> parameters);
		double LogLikelihood(std::vector <double> parameters);

 private:

		BCModel* fTestModel; 

};
// ---------------------------------------------------------

#endif

