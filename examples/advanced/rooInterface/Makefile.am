include $(top_srcdir)/examples/Makefile-common.am

# move outside of ROOSTATS section to avoid
# `examples/Makefile-common.am:21: warning: 'TEST_EXTENSIONS' cannot have conditional contents`
TEST_EXTENSIONS += .root

# Requires BAT build with RooStats enabled
if ROOSTATS

noinst_PROGRAMS = runRooInterface

runRooInterface_SOURCES = runRooInterface.cxx

input_workspaces = \
	WS_GaussOverFlat.root \
	WS_GaussOverFlat_withSystematics.root \
	WS_Poisson.root \
	WS_Poisson_withSystematics.root

CLEANFILES += *.txt
$(input_workspaces): WS_%.root: PrepareWorkspace_%.C
	$(ROOTEXEC) -l -x -q $<

ROOT_LOG_COMPILER = $(builddir)/runRooInterface

TESTS = $(input_workspaces) runBATCalculator.C

endif
