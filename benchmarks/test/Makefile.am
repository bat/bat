CLEANFILES = \
	parallel_TEST_GaussModelSerial.root \
	parallel_TEST_GaussModelParallel.root \
	BCModel_TEST.pdf \
	BCSummaryTool_TEST_summary_par.pdf \
	BCSummaryTool_TEST_summary_corr.pdf \
	BCSummaryTool_TEST_summary_upd.pdf \
	*~ \
	log.txt

MAINTAINERCLEANFILES = Makefile.in

AM_CXXFLAGS = -I$(top_srcdir) $(ROOTCFLAGS)

if DEBUG
AM_CXXFLAGS += -g -O0 -Wall -Wextra -pedantic -Wno-long-long -Wno-unknown-pragmas
else
AM_CXXFLAGS += -O2 -Wno-unknown-pragmas
endif

noinst_LIBRARIES = libBATtest.a

libBATtest_a_SOURCES = \
	test.cxx test.h \
	GaussModel.cxx GaussModel.h \
	IntegrationModel.cxx IntegrationModel.h

TESTS = \
	test_TEST \
	BCMath_TEST \
	BCModel_TEST \
	BCSummaryTool_TEST \
	BCIntegrate_TEST \
	parallel_TEST

LDADD = libBATtest.a $(top_builddir)/src/libBAT.la \
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

BCIntegrate_TEST_SOURCES = BCIntegrate_TEST.cxx

parallel_TEST_SOURCES = parallel_TEST.cxx
