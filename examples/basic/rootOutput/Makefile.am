include $(top_srcdir)/examples/Makefile-common.am

noinst_PROGRAMS = readRootOutput runRootOutput

readRootOutput_SOURCES = readRootOutput.cxx
runRootOutput_SOURCES = runRootOutput.cxx GaussModel.cxx

TESTS = runRootOutput.cxx readRootOutput.cxx

readRootOutput.log: runRootOutput.log
