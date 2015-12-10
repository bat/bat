CLEANFILES = \
	*.pdf \
	*.root \
	*~ \
	log.txt

MAINTAINERCLEANFILES = Makefile.in

AM_CXXFLAGS = $(ROOTAUXCFLAGS) -I$(top_srcdir) $(ROOTCFLAGS)

if DEBUG
AM_CXXFLAGS += -g -O0 -Wall -Wextra -pedantic -Wno-long-long -Wno-unknown-pragmas
else
AM_CXXFLAGS += -O2 -Wno-unknown-pragmas
endif
AM_LDFLAGS = -L$(ROOTLIBDIR)

noinst_LIBRARIES = libBATtest.a

libBATtest_a_SOURCES = \
	test.cxx test.h \
	GaussModel.cxx GaussModel.h

# Use extension .TEST for test binaries need to make test depencencies work:
TEST_EXTENSIONS = .TEST

TESTS = \
	test.TEST \
	BCAux.TEST \
	BCEngineMCMC.TEST \
	BCMath.TEST \
	BCModel.TEST \
	BCParameter.TEST \
	BCPrior.TEST \
	BCSummaryTool.TEST \
	parallel.TEST

# If parallelization is enabled, certain tests that use parallelization
# internally should be run sequentially, not in parallel, to avoid
# massive performance breakdown.
# The order is determined because parallel and BCSummaryTool also run Markov chains
# and we want BCEngineMCMC to test that feature first.
if THREAD_PARALLELIZATION
BCSummaryTool.log : BCEngineMCMC.log
parallel.log : BCSummaryTool.log
endif

LDADD = \
	libBATtest.a \
	$(top_builddir)/src/libBAT.la \
	$(top_builddir)/models/base/libBATmodels.la \
	$(top_builddir)/models/mtf/libBATmtf.la \
	$(ROOTLIBS)

check_PROGRAMS = $(TESTS)

test_TEST_SOURCES = test_TEST.cxx

BCAux_TEST_SOURCES = BCAux_TEST.cxx

BCEngineMCMC_TEST_SOURCES = BCEngineMCMC_TEST.cxx

BCMath_TEST_SOURCES = BCMath_TEST.cxx

BCModel_TEST_SOURCES = BCModel_TEST.cxx

BCParameter_TEST_SOURCES = BCParameter_TEST.cxx

BCPrior_TEST_SOURCES = BCPrior_TEST.cxx

BCSummaryTool_TEST_SOURCES = BCSummaryTool_TEST.cxx

parallel_TEST_SOURCES = parallel_TEST.cxx

# need special Legendre function from the gsl via mathmore
# for the integration test functions
if ROOTMATHMORE
libBATtest_a_SOURCES +=	IntegrationModel.cxx IntegrationModel.h
TESTS += BCIntegrate.TEST
BCIntegrate_TEST_SOURCES = BCIntegrate_TEST.cxx
if THREAD_PARALLELIZATION
BCIntegrate.log: BCSummaryTool.log
endif
endif
