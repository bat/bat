#!/bin/bash

MYPROJECT=$1
MYMODEL=$2

if [[ $2 == '' ]]; then
	MYMODEL=$MYPROJECT
fi

function usage() {
cat << EOF
Usage: $0 <project_name> [class_name]

Creates a new directory <project_name> and a new model class <class_name>,
i.e. the header and source files for the class, together with a Makefile
for the project.

EOF
}

if [[ "$1" == '' ]]; then
	usage
	exit 0
fi

upMYMODEL=`echo $MYMODEL|tr [a-z] [A-Z]`

function m_header() {
cat << EOF
// ***************************************************************
// This file was created using the $0 script
// for project $MYPROJECT
// $0 is part of Bayesian Analysis Toolkit (BAT).
// BAT can be downloaded from http://www.mppmu.mpg.de/bat
// ***************************************************************

#ifndef __|:UP_MODEL:|__H
#define __|:UP_MODEL:|__H

#include <BAT/BCModel.h>

// This is a |:Model:| header file.
// Model source code is located in file |:Project:|/|:Model:|.cxx

// ---------------------------------------------------------
class |:Model:| : public BCModel
{
	public:

		// Constructors and destructor
		|:Model:|();
		|:Model:|(const char * name);
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
// ***************************************************************
// This file was created using the $0 script
// for project $MYPROJECT
// $0 is part of Bayesian Analysis Toolkit (BAT).
// BAT can be downloaded from http://www.mppmu.mpg.de/bat
// ***************************************************************

#include "|:Model:|.h"

// ---------------------------------------------------------
|:Model:|::|:Model:|() : BCModel()
{  // default constructor
	DefineParameters();
};

// ---------------------------------------------------------
|:Model:|::|:Model:|(const char * name) : BCModel(name)
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
	// You can then use them in the methods below by calling the
	// parameters.at(i) or parameters[i], where i is the index
	// of the parameter. The indices increase from 0 according to the
	// order of adding the parameters.

//	this -> AddParameter("par1", 0.0, 1.0);   // index 0
//	this -> AddParameter("par2", -7.0, 28.0); // index 1
}

// ---------------------------------------------------------
double |:Model:|::LogLikelihood(std::vector <double> parameters)
{
	// This methods returns the logarithm of the conditional probability
	// p(data|parameters). This is where you have to define your model.

	double logprob = 0.;

	return logprob;
}

// ---------------------------------------------------------
double |:Model:|::LogAPrioriProbability(std::vector <double> parameters)
{
	// This method returns the logarithm of the prior probability for the
	// parameters p(parameters).

	double logprob = 0.;

	// For flat prior it's very easy.
//	for(unsigned int i=0; i < this -> GetNParameters(); i++)
//		logprob -= log(this -> GetParameter(i) -> GetRangeWidth());

	return logprob;
}
// ---------------------------------------------------------

EOF
}

function m_makefile() {
cat << EOF
###################################################################
# This Makefile was created using the $0 script
# for project $MYPROJECT
# $0 is part of Bayesian Analysis Toolkit (BAT).
# BAT can be downloaded from http://www.mppmu.mpg.de/bat
###################################################################
#
# Run 'make' to compile the program and 'make clean' to remove
# all compiled parts and 'clean' the directory.
#
# You might need to adjust the CFLAGS, LIBS, and GLIBS based on
# the BAT installation on your system. Consult the gmake manual
# for details.
#
###################################################################

# Root variables
ROOTCFLAGS   := \$(shell root-config --cflags)
ROOTLIBS     := \$(shell root-config --libs) -lMinuit
ROOTGLIBS    := \$(shell root-config --glibs)

# compiler and flags
CXX          = g++
CXXFLAGS     = -g -Wall -fPIC -Wno-deprecated -O2
LD           = g++
LDFLAGS      = -g -O2
SOFLAGS      = -shared

# standard commands
RM           = rm -f
MV           = mv
ECHO         = echo
CINT         = rootcint

# add ROOT flags
CXXFLAGS    += \$(ROOTCFLAGS)

# ----------------------------------------------------------------------
# The following definitions depend on the setup of the system where
# the project is being compiled. If BAT is installed in the standard
# system search path then the lines below are correct and the compilation
# will work
CXXFLAGS    += -I. -I./include
LIBS        += \$(ROOTLIBS)  -lBAT
GLIBS       += \$(ROOTGLIBS) -lBAT

# In case you see following errors during the compilation,
#
#   undefined reference to 'Divonne'
#   undefined reference to 'Suave'
#   undefined reference to 'Cuhre'
#   undefined reference to 'Vegas'
#
# your version of BAT was installed with the Cuba support and you need
# to adjust the Makefile by uncommenting the following lines. You might
# also need to add the path to libcuba.a to the lines below
# as -L/path/to/cuba/lib
#
# LIBS        += -lcuba
# GLIBS       += -lcuba

# If BAT was installed in a non-standard path (e.g. user home
# directory) then one can specify the location of the header
# files here (uncomment following lines and adjust the path in
# BATINSTALLDIR):
#
# BATINSTALLDIR = '/path/to/bat/installation/directory'
# CXXFLAGS    += -I\$(BATINSTALLDIR)/include
# LIBS        += -L\$(BATINSTALLDIR)/lib -lBAT
# GLIBS       += -L\$(BATINSTALLDIR)/lib -lBAT

# List of all classes (models) used in the program
# Add classes to the end. Baskslash indicates continuation
# on the next line
CXXSRCS      = \\
        |:Model:|.cxx

# ----------------------------------------------------------------------
# don't change lines below unless you know what you're doing
#

CXXOBJS      = \$(patsubst %.cxx,%.o,\$(CXXSRCS))
EXEOBJS      =
MYPROGS     = \\
        run|:Project:|

GARBAGE      = \$(CXXOBJS) \$(EXEOBJS) *.o *~ link.d \$(MYPROGS)


# targets
all : project

link.d : \$(patsubst %.cxx,%.h,\$(CXXSRCS))
	\$(CXX) -MM \$(CXXFLAGS) \$(CXXSRCS) > link.d;

include link.d

%.o : %.cxx
	\$(CXX) \$(CXXFLAGS) -c \$< -o \$@

clean :
	\$(RM) \$(GARBAGE)

project : run|:Project:|.cxx \$(CXXOBJS)
	\$(CXX) \$(CXXFLAGS) -c \$<
	\$(CXX) \$(LDFLAGS) \$(LIBS) run|:Project:|.o \$(CXXOBJS) -o run|:Project:|

print :
	echo compiler  : \$(CXX)
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
// ***************************************************************
// This file was created using the $0 script
// for project $MYPROJECT
// $0 is part of Bayesian Analysis Toolkit (BAT).
// BAT can be downloaded from http://www.mppmu.mpg.de/bat
// ***************************************************************

#include <BAT/BCLog.h>
#include <BAT/BCAux.h>

#include "|:Model:|.h"

int main()
{

	// set nice style for drawing than the ROOT default
	BCAux::SetStyle();

	// open log file
//	BCLog::OpenLog("log.txt");
	BCLog::SetLogLevel(BCLog::detail);

	// create new |:Model:| object
	|:Model:| * m = new |:Model:|();

	BCLog::OutSummary("Test model created");

	// perform your analysis here

	// normalize the posterior, i.e. integrate posterior
	// over the full parameter space
//	m -> Normalize();

	// run MCMC and marginalize posterior wrt. all parameters
	// and all combinations of two parameters
//	m -> MarginalizeAll();

	// run mode finding; by default using Minuit
//	m -> FindMode();

	// if MCMC was run before (MarginalizeAll()) it is
	// possible to use the mode found by MCMC as
	// starting point of Minuit minimization
//	m -> FindMode( m -> GetBestFitParameters() );

	// draw all marginalized distributions into a PostScript file
//	m -> PrintAllMarginalized("|:Model:|_plots.ps");

	// calculate p-value
//	m -> CalculatePValue( m -> GetBestFitParameters() );

	// print results of the analysis into a text file
//	m -> PrintResults("|:Model:|_results.txt");

	// close log file
//	BCLog::CloseLog();

	delete m;

	BCLog::OutSummary("Test program ran successfully");
	BCLog::OutSummary("Exiting");

	return 0;

}

EOF
}

mkdir $MYPROJECT

m_header | sed -e "s/|\:UP_MODEL\:|/$upMYMODEL/g; s/|\:Model\:|/$MYMODEL/g; s/|\:Project\:|/$MYPROJECT/g" > $MYPROJECT/$MYMODEL.h
m_sourcecode | sed -e "s/|\:Model\:|/$MYMODEL/g" > $MYPROJECT/$MYMODEL.cxx
m_makefile | sed -e "s/|\:Model\:|/$MYMODEL/g; s/|\:Project\:|/$MYPROJECT/g" > $MYPROJECT/Makefile
m_main | sed -e "s/|\:Model\:|/$MYMODEL/g; s/|\:Project\:|/$MYPROJECT/g" > $MYPROJECT/run$MYPROJECT.cxx

cat << EOF
--------------------------------------------------------------------------
The new BAT project was created in the directory '$MYPROJECT'.
To test the configuration try to compile the project by running 'make'
inside the directory. In case there are some compilation errors you need
to adjust the parameters inside the 'Makefile'.

Once the program is compiled successfully, you can run it and it should
print some basic information on the screen.

Implement your model in file:    $MYMODEL.cxx
Implement your analysis in file: $MYPROJECT.cxx

Consult BAT webpage for details: http://www.mppmu.mpg.de/bat
--------------------------------------------------------------------------

EOF
