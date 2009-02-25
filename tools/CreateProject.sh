#!/bin/bash

MYPROJECT=$1
MYMODEL=$2

if [[ $2 == '' ]]; then
	MYMODEL=$MYPROJECT
fi

function usage() {
cat << EOF
Usage: $0 <project_name> [class_name]

Script create a new directory <project_name> and creates a new model class,
i.e. the header and source files for the class.

EOF
}

if [[ "$1" == '' ]]; then
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
// Model source code is located in file |:Project:|/src/|:Model:|.cxx

// ---------------------------------------------------------
class |:Model:| : public BCModel
{
	public:

		// Constructors and destructor
		|:Model:|();
		|:Model:|(const char* name);
		~|:Model:|();

		// Methods to overload, see file |:Model:|.cxx
		void DefineParameters();
		double LogAPrioriProbability(std::vector <double> parameters);
		double LogLikelihood(std::vector <double> parameters);
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
	// add parameters to your model here.
	//
	// you can then use them in the methods below by calling the
	// parameters.at(i) or parameters[i], where i is the index
	// of the parameter. the indices increase from 0 according to the
	// order of adding the parameters.

//	this -> AddParameter("par1", 0.0, 1.0);   // index 0
//	this -> AddParameter("par2", -7.0, 28.0); // index 1
}

// ---------------------------------------------------------
double |:Model:|::LogLikelihood(std::vector <double> parameters)
{
	// this methods returns the logarithm of the conditional probability
	// p(data|parameters).

	double logprob = 0.;

	return logprob;
}

// ---------------------------------------------------------
double |:Model:|::LogAPrioriProbability(std::vector <double> parameters)
{

	// this method returs the logarithm of the prior probability for the
	// parameters p(parameters).

	double logprob = 0.;

// logprob -= fabs(0.5 - parameters.at(0));
// logprob -= 0.1*(parameters.at(1)-10.)*(parameters.at(1)-10.);

	return logprob;
}
// ---------------------------------------------------------

EOF
}

function m_makefile() {
cat << EOF
# Root variables
ROOTCFLAGS   := \$(shell root-config --cflags)
ROOTLIBS     := \$(shell root-config --libs) -lMinuit
ROOTGLIBS    := \$(shell root-config --glibs)

# Programs
CXX          = g++
CXXFLAGS     = -g -Wall -fPIC -Wno-deprecated -O2
LD           = g++
LDFLAGS      = -g -O2
SOFLAGS      = -shared

RM           = rm -f
MV           = mv
ECHO         = echo
CINT         = rootcint

# Assign or Add variables
CXXFLAGS    += \$(ROOTCFLAGS)
CXXFLAGS    += -I. -I./include -I\$(BATINSTALLDIR)/BAT

LIBS        += \$(ROOTLIBS)  -L\$(BATINSTALLDIR)/lib -lBAT
GLIBS       += \$(ROOTGLIBS) -L\$(BATINSTALLDIR)/lib -lBAT

CXSRCS      = \
        |:Model:|.cxx

CXXSRCS      = \$(patsubst %.cxx,src/%.cxx,\$(CXSRCS))
CXXOBJS      = \$(patsubst %.cxx,%.o,\$(CXSRCS))
EXEOBJS      =
GARBAGE      = \$(CXXOBJS) \$(EXEOBJS) *.o *~ link.d |:Project:|.exe

all : project

link.d : \$(patsubst %.cxx,include/%.h,\$(CXSRCS))
	\$(CXX) -MM \$(CXXFLAGS) \$(CXXSRCS) > link.d;

include link.d

%.o : src/%.cxx
	\$(CXX) \$(CXXFLAGS) -c \$< -o \$@

clean :
	\$(RM) \$(GARBAGE)

project : run|:Project:|.cxx \$(CXXOBJS)
	\$(CXX) \$(CXXFLAGS) -c \$<
	\$(CXX) \$(LDFLAGS) \$(LIBS) run|:Project:|.o \\
		\$(CXXOBJS) -o run|:Project:|.exe

print :
	echo compiler  : \$(CXX)
	echo compiler  : \$(CXSRCS)
	echo c++ srcs  : \$(CXXSRCS)
	echo c++ objs  : \$(CXXOBJS)
	echo c++ flags : \$(CXXFLAGS)
	echo libs      : \$(LIBS)
	echo so flags  : \$(SOFLAGS)

	echo rootlibs  : \$(ROOTLIBS)
	echo rootglibs : \$(ROOTGLIBS)

EOF
}

function m_main() {
cat << EOF
#include "|:Model:|.h"
#include <BAT/BCLog.h>
#include <BAT/BCAux.h>

int main()
{

	// set nice style for drawing than the ROOT default
	BCAux::SetStyle();

	// open log file with default level of logging
	BCLog::OpenLog("log.txt");
	BCLog::SetLogLevel(BCLog::detail);

	// create new |:Model:| object
	|:Model:| * _|:Model:| = new |:Model:|();

	// perform your analysis here
//	_|:Model:| -> Normalize();
//	_|:Model:| -> MarginalizeAll();
//	_|:Model:| -> FindMode( _|:Model:| -> GetBestFitParameters() );
//	_|:Model:| -> PrintAllMarginalized("plots.ps");
//	_|:Model:| -> CalculatePValue( _|:Model:| -> GetBestFitParameters() );
//	_|:Model:| -> PrintResults("results.txt");

	// close log file
	BCLog::CloseLog();

	return 0;

}

EOF
}

echo "Create project" $MYPROJECT
mkdir $MYPROJECT
mkdir $MYPROJECT/src
mkdir $MYPROJECT/include

m_header | sed -e "s/|\:UP_MODEL\:|/$upMYMODEL/g; s/|\:Model\:|/$MYMODEL/; s/|\:Project\:|/$MYPROJECT/" > $MYPROJECT/include/$MYMODEL.h
m_sourcecode | sed -e "s/|\:Model\:|/$MYMODEL/g" > $MYPROJECT/src/$MYMODEL.cxx
m_makefile | sed -e "s/|\:Model\:|/$MYMODEL/g; s/|\:Project\:|/$MYPROJECT/" > $MYPROJECT/Makefile
m_main | sed -e "s/|\:Model\:|/$MYMODEL/g; s/|\:Project\:|/$MYPROJECT/" > $MYPROJECT/run$MYPROJECT.cxx
