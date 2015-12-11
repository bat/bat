include $(top_srcdir)/examples/Makefile-common.am

# move outside of ROOSTATS section to avoid
# `examples/Makefile-common.am:21: warning: 'TEST_EXTENSIONS' cannot have conditional contents`
TEST_EXTENSIONS += .root

# Requires BAT build with RooStats enabled
if ROOSTATS
# examples are not thread safe
if THREAD_PARALLELIZATION
AM_TESTS_ENVIRONMENT += export OMP_NUM_THREADS=1;
endif

noinst_PROGRAMS = runRooInterface

runRooInterface_SOURCES = runRooInterface.cxx

input_workspaces = \
	WS_GaussOverFlat.root \
	WS_GaussOverFlat_withSystematics.root \
	WS_Poisson.root \
	WS_Poisson_withSystematics.root

CLEANFILES += *.txt
$(input_workspaces): WS_%.root: PrepareWorkspace_%.C
	$(C_LOG_COMPILER) $(AM_C_LOG_FLAGS) $<

ROOT_LOG_COMPILER = $(builddir)/runRooInterface

TESTS = $(input_workspaces)

# import in macro only works with root 5
if !WITH_CLING
TESTS += runBATCalculator.C
endif

endif
