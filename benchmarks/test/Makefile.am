CLEANFILES = \
	*~ \
	log.txt
	
MAINTAINERCLEANFILES = Makefile.in

AM_CXXFLAGS = -I$(top_srcdir)

if DEBUG
AM_CXXFLAGS += -g -O0 -Wextra -pedantic -Wno-long-long
else !DEBUG
AM_CXXFLAGS += -O2
endif !DEBUG

noinst_LIBRARIES = libBATtest.a

libBATtest_a_SOURCES = \
	test.cxx test.h

TESTS = \
	test_TEST \
	BCMath_TEST

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
