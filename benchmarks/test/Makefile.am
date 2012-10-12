CLEANFILES = \
	parallel_TEST_GaussModelSerial.root \
	parallel_TEST_GaussModelParallel.root \
	*~ \
	log.txt

MAINTAINERCLEANFILES = Makefile.in

AM_CXXFLAGS = -std=c++0x -I$(top_srcdir) $(ROOTCFLAGS)

if DEBUG
AM_CXXFLAGS += -g -O0 -Wall -Wextra -pedantic -Wno-long-long -Wno-unknown-pragmas
else !DEBUG
AM_CXXFLAGS += -O2 -Wno-unknown-pragmas
endif !DEBUG

if THREAD_PARALLELIZATION
AM_CXXFLAGS += -fopenmp
endif

noinst_LIBRARIES = libBATtest.a

libBATtest_a_SOURCES = \
	test.cxx test.h \
	GaussModel.cxx GaussModel.h

TESTS = \
	test_TEST \
	BCMath_TEST \
	parallel_TEST

TESTS_LD_ADD = libBATtest.a $(top_builddir)/src/libBAT.la $(ROOTLIBS) -lMinuit
TESTS_LD_FLAGS = -L$(ROOTLIBDIR)

if CUBA
TESTS_LD_ADD += -lcuba
endif

check_PROGRAMS = $(TESTS)

test_TEST_SOURCES = test_TEST.cxx
test_TEST_LDADD = $(TESTS_LD_ADD)
test_TEST_LDFLAGS = $(TESTS_LD_FLAGS)

BCMath_TEST_SOURCES = BCMath_TEST.cxx
BCMath_TEST_LDADD = $(TESTS_LD_ADD)
BCMath_TEST_LDFLAGS = $(TESTS_LD_FLAGS)

parallel_TEST_SOURCES = parallel_TEST.cxx
parallel_TEST_LDADD = $(TESTS_LD_ADD)
parallel_TEST_LDFLAGS = $(TESTS_LD_FLAGS)