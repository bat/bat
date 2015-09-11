include $(top_srcdir)/examples/Makefile-common.am

CLEANFILES += data.txt

noinst_PROGRAMS = asymmErrorsFitterExample

asymmErrorsFitterExample_SOURCES = asymmErrorsFitterExample.cxx PolAsymm.cxx

TESTS = CreateData.C asymmErrorsFitterExample.cxx

asymmErrorsFitterExample.log: CreateData.log
