#include "BCParameter.h" 

// --------------------------------------------------------- 

BCParameter::BCParameter()
{

	fName       = "parameter"; 
	fLowerLimit = 0.0; 
	fUpperLimit = 1.0; 
	fNuisance   = 0;

}

// --------------------------------------------------------- 

BCParameter::BCParameter(const char * name, double lowerlimit, double upperlimit)
{

	BCParameter(); 

	fName       = name; 
	fLowerLimit = lowerlimit; 
	fUpperLimit = upperlimit; 

}

// --------------------------------------------------------- 

BCParameter::BCParameter(const BCParameter & parameter)
{

	parameter.Copy(*this); 

}

// --------------------------------------------------------- 

BCParameter & BCParameter::operator = (const BCParameter & parameter)
{

	if (this != &parameter) 
		parameter.Copy(* this); 

	return * this; 

}

// --------------------------------------------------------- 

BCParameter::~BCParameter()
{

}

// --------------------------------------------------------- 

void BCParameter::PrintSummary()
{

	cout
		<<"       > Parameter   : "<< fName <<endl
		<<"         Index       : "<< fIndex <<endl
		<<"         Lower Limit : "<< fLowerLimit <<endl
		<<"         Upper Limit : "<< fUpperLimit <<endl
		<<"         Nuisance    : ";

	if(fNuisance)
		cout<<"Yes"<<endl; 
	else
		cout<<"No"<<endl; 

	cout << endl;

}

// --------------------------------------------------------- 

void BCParameter::Copy(BCParameter & parameter) const 
{

	parameter.fName       = this -> fName; 
	parameter.fIndex      = this -> fIndex; 
	parameter.fLowerLimit = this -> fLowerLimit;  
	parameter.fUpperLimit = this -> fUpperLimit; 
	parameter.fNuisance   = this -> fNuisance;

}

// --------------------------------------------------------- 
