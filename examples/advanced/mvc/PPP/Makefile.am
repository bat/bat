include $(top_srcdir)/examples/Makefile-common.am

CLEANFILES += MyCombination_full_BLUE.txt

noinst_PROGRAMS = runMyCombination

runMyCombination_SOURCES = runMyCombination.cxx MyCombination.cxx

TEST_EXTENSIONS += .txt
TXT_LOG_COMPILER = $(builddir)/runMyCombination

# Disabled:
# TESTS = input.txt
