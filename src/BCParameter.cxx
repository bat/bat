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

	cout
		<<"       > Parameter   : "<< fName <<endl
		<<"         Index       : "<< fIndex <<endl
		<<"         Lower Limit : "<< fLowerLimit <<endl
		<<"         Upper Limit : "<< fUpperLimit <<endl
		<<"         Nuisence    : ";

	if(fNuisence)
		cout<<"Yes"<<endl; 
	else
		cout<<"No"<<endl; 

	cout << endl;

}

// --------------------------------------------------------- 
