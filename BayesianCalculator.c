#include "BCModel1.h" 
#include "BCModel2.h" 
#include "BCModel.h"
#include "BCModelParameter.h" 
#include "BCManager.h" 
#include "BCDataObject.h" 

#include <TH1F.h>
#include  <iostream.h>

// ---------------------------------------------------------
  
int main()
{

  // ---------------------------------------------------------
  // define model1 
  // ---------------------------------------------------------

  BCModel1* fModel1 = new BCModel1("PeakModel"); 

  // define parameters 

  BCModelParameter* model1_parameter1 = new BCModelParameter("events", 0.0, 40.0); 
  BCModelParameter* model1_parameter2 = new BCModelParameter("centralmass", 170.0, 180.0); 

  // add parameters to model 

  fModel1 -> AddParameter(model1_parameter1); 
  fModel1 -> AddParameter(model1_parameter2); 

  // ---------------------------------------------------------
  // define model2 
  // ---------------------------------------------------------

  BCModel2* fModel2 = new BCModel2("FlatModel"); 

  // define parameters 

  BCModelParameter* model2_parameter1 = new BCModelParameter("events", 0.0, 40.0); 
  
  // add parameters to model 

  fModel2 -> AddParameter(model2_parameter1); 

  // ---------------------------------------------------------
  // define manager 
  // ---------------------------------------------------------

  BCManager* fManager = new BCManager(); 

  // ---------------------------------------------------------
  // add models to model manager with a priori probabilities 
  // ---------------------------------------------------------

  fManager -> AddModel(fModel1, 0.2); 
  fManager -> AddModel(fModel2, 0.8); 

  // ---------------------------------------------------------
  // read data from file 
  // ---------------------------------------------------------

  // read data from ROOT file 
  
  fManager -> ReadTreeFromRootFile("tools/test.root", "testtree", "mass,energy"); 
  
  // ---------------------------------------------------------
  // initialize 
  // ---------------------------------------------------------

  fManager -> Initialize(); 

  // ---------------------------------------------------------
  // summarize
  // ---------------------------------------------------------
  
  fManager -> PrintSummary(); 

  // ---------------------------------------------------------
  // Marginalize probabilities 
  // ---------------------------------------------------------

  fModel1 -> MarginalizeProbabilities(); 

  return 0; 

}

// ---------------------------------------------------------
  
