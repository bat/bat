include $(top_srcdir)/examples/Makefile-common.am

CLEANFILES += systematics/*.pdf

noinst_PROGRAMS = runSystematics

runSystematics_SOURCES = runSystematics.cxx

TESTS = CreateHistograms.C runSystematics.cxx

runSystematics.log: CreateHistograms.log
