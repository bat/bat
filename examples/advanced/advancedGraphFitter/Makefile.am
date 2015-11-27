include $(top_srcdir)/examples/Makefile-common.am

noinst_PROGRAMS =

TESTS = graphFitterAdvancedExample.C

# the input data is only in $(srcdir) but for an out-of-source build
# we need access. `ln -s` fails for in-source build because data.txt
# already exists and that's fine.
$(abs_builddir)/data.txt :
	ln -s $(srcdir)/data.txt data.txt

graphFitterAdvancedExample.log : $(abs_builddir)/data.txt
