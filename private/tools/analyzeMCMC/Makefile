ROOTCFLAGS   := $(shell root-config --cflags)
ROOTLIBS     := $(shell root-config --libs) -lMinuit
ROOTGLIBS    := $(shell root-config --glibs)

CXX          = g++
CXXFLAGS     = -g -Wall -fPIC -Wno-deprecated -O2 $(ROOTCFLAGS)
all: 
	$(CXX) $(CXXFLAGS) $(ROOTLIBS) $(ROOTGLIBS) analyzeMCMC.c -o analyzeMCMC.exe
