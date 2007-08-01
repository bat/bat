#!/bin/bash

MYMODEL=$1

function usage() {
cat << EOF
Usage: $0 <class_name> 
Script takes one argument which is a name of the new model class
and creates header and source files for the class.

EOF
}

if [[ "$1"x == x ]]; then
	usage
	exit 0
fi

upMYMODEL=`echo $MYMODEL|tr [a-z] [A-Z]`

function m_header() {
cat << EOF
#ifndef __|:UP_MODEL:|__H
#define __|:UP_MODEL:|__H

#include "BCModel.h" 

// This is a |:Model:| header file.
// Model source code is located in file |:Model:|.cxx

// --------------------------------------------------------- 
class |:Model:| : public BCModel 
{
  public: 

  // Coonstructors and destructor
    |:Model:|();
    |:Model:|(const char* name); 
    ~|:Model:|();


  // Methods to overload, see file |:Model:|.cxx

    void DefineParameters(); 

    double APrioriProbability(std::vector <double> parameters); 

    double Likelihood(std::vector <double> parameters);

//    double ConditionalProbabilityEntry(BCDataPoint* datapoint, std::vector <double> parameters); 

//    double PoissonProbability(int nentries, std::vector <double> parameters); 

}; 
// --------------------------------------------------------- 

#endif 

EOF
}

function m_sourcecode() {
cat << EOF
#include "|:Model:|.h" 

// --------------------------------------------------------- 
|:Model:|::|:Model:|() : BCModel()
{  // default constructor
   DefineParameters();
};

// --------------------------------------------------------- 
|:Model:|::|:Model:|(const char* name) : BCModel(name)
{  // constructor
   DefineParameters();
};

// --------------------------------------------------------- 
|:Model:|::~|:Model:|()
{};  // default destructor

// --------------------------------------------------------- 
void |:Model:|::DefineParameters()
{
// Add parameters to your model here.
//
// You can then use them in the methods below by calling the
// parameters.at(i) or parameters[i], where i is the index
// of the parameter. The indeces increase from 0 according to the
// order of adding the parameters.

//   this -> AddParameter("par1", 0.0, 1.0);   // index 0 
//   this -> AddParameter("par2", -7.0, 28.0); // index 1

}

// --------------------------------------------------------- 
double |:Model:|::Likelihood(std::vector <double> parameters)
{
// Put the code for calculating the Likelyhood p(data|params) here.
// Alternatively you can overload the methods ConditionalProbabilityEntry
// and PoissonProbability below and comment out this method.
// See BC documentation for more details.

// This method is called many times during integration and marginalization,                      
// so it should be rather efficient. I.e. it should use parameters[i]                            
// instead of calls to GetParameter("name"). The user should make sure that 
// he calls the parameters by the index correctly.                                               


// The likelihood does not have to be normalized.

   double probability = 1.0; 

   return probability; 
}

// --------------------------------------------------------- 
double |:Model:|::APrioriProbability(std::vector <double> parameters)
{
// Put the code for calculating the Prior probability p_0(params) here.
// The prior does not have to be normalized.

   double probability = 1.0; 

   return probability; 
}

// ---------------------------------------------------------                                     
// ---------------------------------------------------------                                     
// Alternatively to overloading the Likelihood method, you can overload                          
// the two methods below. In this case you have to comment out the                               
// Likelihood method in this file as well as in the header file.
// If you use these methods the likelihood is then calculated as a product                       
// of ConditionalProbabilityEntry() for individual data points times
// PoissonProbability(). See BC documentation for more details.    
//
/*                                                                                               
// --------------------------------------------------------- 
double |:Model:|::ConditionalProbabilityEntry(BCDataPoint* datapoint,
          std::vector <double> parameters)
{
   double probability = 1.0; 

   return probability; 
}

// --------------------------------------------------------- 
double |:Model:|::PoissonProbability(int nentries,
          std::vector <double> parameters)
{
   double probability = 1.0; 

   return probability; 
}
*/

// --------------------------------------------------------- 

EOF
}

m_header | sed -e "s/|\:UP_MODEL\:|/$upMYMODEL/g; s/|\:Model\:|/$MYMODEL/" > $MYMODEL.h
m_sourcecode | sed -e "s/|\:Model\:|/$MYMODEL/g" > $MYMODEL.cxx

