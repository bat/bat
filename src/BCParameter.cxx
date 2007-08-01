#include "BCParameter.h" 

// --------------------------------------------------------- 

BCParameter::BCParameter()
{

  fName = new char; 
  fName = "parameter"; 
  fLowerLimit = 0.0; 
  fUpperLimit = 1.0; 
  fNuisence = 0;
}

// --------------------------------------------------------- 

BCParameter::BCParameter(const char* name, double lowerlimit, double upperlimit)
{

  fName = (char*) name; 
  fLowerLimit = lowerlimit; 
  fUpperLimit = upperlimit; 

}

// --------------------------------------------------------- 

BCParameter::~BCParameter()
{

  delete fName; 

}

// --------------------------------------------------------- 

void BCParameter::PrintSummary()
{

  std::cout << "     Parameter                  : " << fName << std::endl; 
  std::cout << "     Index                      : " << fIndex << std::endl; 
  std::cout << "     Lower Limit                : " << fLowerLimit << std::endl; 
  std::cout << "     Upper Limit                : " << fUpperLimit << std::endl; 
  if(fNuisence)
	  std::cout << "     Is Nuisence Parameter."<< std::endl; 
  else
	  std::cout << "     Is Not Nuisence Parameter."<< std::endl; 
  std::cout << std::endl;

}

// --------------------------------------------------------- 


