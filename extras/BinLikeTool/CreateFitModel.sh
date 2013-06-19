#!/bin/bash

MYPROJECT=$1
MYMODEL=$2

if [[ $2 == '' ]]; then
	MYMODEL=$MYPROJECT
fi

function usage() {
cat << EOF
Usage: $0 [class_name]

Creates a new model class <class_name>, i.e. the header and source
files for the class, and a run file. 

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
// $0 is part of Bayesian Analysis Toolkit (BAT).
// BAT can be downloaded from http://mpp.mpg.de/bat
// ***************************************************************

#ifndef __|:UP_MODEL:|__H
#define __|:UP_MODEL:|__H

#include <../BinLikeModel.h>

// This is a |:Model:| header file.
// Model source code is located in file |:Project:|/|:Model:|.cxx

// ---------------------------------------------------------
class |:Model:| : public BinLikeModel
{
  public:

  /**
   * The default constructor.
   */
  |:Model:|();

  /**
   * The default destructor.
   */
  ~|:Model:|();

  /**
   * Defines the parameters of the fit. */ 
  void DefineParameters(); 

  /**
   * Calculates the expectation value. 
   * @param parameters The current parameter values.
   * @param parvalue the leading parameter value. 
   * @param x the x-value
   * @return The expectation value. */ 
  double Expectation(std::vector <double> parameters, double parvalue, double x); 
};
// ---------------------------------------------------------

#endif

EOF
}

function m_sourcecode() {
cat << EOF
// ***************************************************************
// This file was created using the $0 script
// $0 is part of Bayesian Analysis Toolkit (BAT).
// BAT can be downloaded from http://mpp.mpg.de/bat
// ***************************************************************

#include <|:Model:|.h>
#include <TMath.h>
// ---------------------------------------------------------
|:Model:|::|:Model:|() : BinLikeModel()
{  
  DefineParameters();
};

// ---------------------------------------------------------
|:Model:|::~|:Model:|()
{}; 

// ---------------------------------------------------------
void |:Model:|::DefineParameters()
{
  // add model parameters 
  //  this -> AddParameter("par1", 0.0, 100.0);   // index 0
  //  this -> AddParameter("par2", 0.0, 100.0); // index 1
}

// ---------------------------------------------------------
double |:Model:|::Expectation(std::vector <double> parameters, double parvalue, double x)
{
  // initialize expectation
  double expectation = 0; 

  // calculate parameters of fit function from model parameters and
  // leading parameter
  // ...  

  // return expectation 
  return expectation;
}
// ---------------------------------------------------------

EOF
}

function m_makefile() {
cat << EOF
###################################################################
# This Makefile was created using the $0 script
# $0 is part of Bayesian Analysis Toolkit (BAT).
# BAT can be downloaded from http://mpp.mpg.de/bat
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
CXXFLAGS    += -I. -I../
LIBS        += \$(ROOTLIBS)  -lBAT -L.. -lBinLikeModel
GLIBS       += \$(ROOTGLIBS) -lBAT -L.. -lBinLikeModel

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

-include link.d

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
// $0 is part of Bayesian Analysis Toolkit (BAT).
// BAT can be downloaded from http://mpp.mpg.de/bat
// ***************************************************************

#include <BAT/BCLog.h>
#include <BAT/BCAux.h>

#include <|:Model:|.h>

#include <TFile.h>
#include <TH1D.h>

int main()
{
	// ----------------------------------------------------
	// configure BAT
	// ----------------------------------------------------
	
	// set nice style for drawing than the ROOT default
	BCAux::SetStyle();

	// open log file
	BCLog::OpenLog("log.txt");
	BCLog::SetLogLevel(BCLog::detail);
	
	// ----------------------------------------------------
	// read histograms from file
	// ----------------------------------------------------

  //	TFile * file = new TFile("input.root", "read"); 
	
  //	TH1D* hist_0 = (TH1D*) file -> Get("hist_0");

	// ----------------------------------------------------
	// create new fit model and perform the fit 
	// ----------------------------------------------------

	// create new fit model 
	|:Model:| * m = new |:Model:|();

	// add histograms
  //	m -> AddHistogram(hist_0, 150.0); 

	// perform fit 
	m -> MarginalizeAll(); 
	m -> FindMode(); 

	// print results
	m -> PrintAllMarginalized("plots.ps"); 
	m -> PrintResults("results.txt"); 
	m -> PrintHistograms(); 
	m -> PrintSummary(); 

	// print goodness-of-fit
	std::cout << " chi2 / dof (chi2-prob) : " 
						<< m->CalculateChi2() << "/" 
						<< m->GetNDF() << " ("
						<< m->CalculateChi2Prob() << ")" << std::endl; 

	// ----------------------------------------------------
	// clean up 
	// ----------------------------------------------------

	// delete model
	delete m; 

	// close file 
  //	file->Close();
  //	delete file; 

	// no errors 
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

Consult BAT webpage for details: http://mpp.mpg.de/bat
--------------------------------------------------------------------------

EOF
