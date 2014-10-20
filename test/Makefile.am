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

noinst_LIBRARIES = libBATtest.a

libBATtest_a_SOURCES = \
	test.cxx test.h \
	GaussModel.cxx GaussModel.h

TESTS = \
	test_TEST \
	BCMath_TEST \
	BCModel_TEST \
	BCSummaryTool_TEST \
	parallel_TEST

LDADD = \
	libBATtest.a \
	$(top_builddir)/src/libBAT.la \
	$(top_builddir)/models/base/libBATmodels.la \
	$(top_builddir)/models/mvc/libBATmvc.la \
	$(top_builddir)/models/mtf/libBATmtf.la \
	$(ROOTLIBS)
	
LDFLAGS = -L$(ROOTLIBDIR)

if CUBA
LDADD += -lcuba
LDFLAGS += @CUBALIB@
endif

check_PROGRAMS = $(TESTS)

test_TEST_SOURCES = test_TEST.cxx

BCMath_TEST_SOURCES = BCMath_TEST.cxx

BCModel_TEST_SOURCES = BCModel_TEST.cxx

BCSummaryTool_TEST_SOURCES = BCSummaryTool_TEST.cxx

parallel_TEST_SOURCES = parallel_TEST.cxx

# need special Legendre function from the gsl via mathmore
# for the integration test functions
if ROOTMATHMORE
libBATtest_a_SOURCES +=	IntegrationModel.cxx IntegrationModel.h
TESTS += BCIntegrate_TEST
BCIntegrate_TEST_SOURCES = BCIntegrate_TEST.cxx
LDADD += -lMathMore
endif
