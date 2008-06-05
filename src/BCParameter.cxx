
#include <iostream>
#include <fstream>

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

	std::cout
		<<"       > Parameter   : "<< fName << std::endl
		<<"         Index       : "<< fIndex << std::endl
		<<"         Lower Limit : "<< fLowerLimit << std::endl
		<<"         Upper Limit : "<< fUpperLimit << std::endl
		<<"         Nuisance    : ";

	if(fNuisance)
		std::cout<<"Yes"<<std::endl; 
	else
		std::cout<<"No"<<std::endl; 

	std::cout << std::endl;

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
