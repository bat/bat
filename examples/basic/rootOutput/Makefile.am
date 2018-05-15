include $(top_srcdir)/examples/Makefile-common.am

noinst_PROGRAMS = readRootOutput runRootOutput continueMCMC

readRootOutput_SOURCES = readRootOutput.cxx
runRootOutput_SOURCES = runRootOutput.cxx GaussModel.cxx
continueMCMC_SOURCES = continueMCMC.cxx GaussModel.cxx

TESTS = runRootOutput.cxx readRootOutput.cxx continueMCMC.cxx

readRootOutput.log: runRootOutput.log
continueMCMC.log: runRootOutput.log
