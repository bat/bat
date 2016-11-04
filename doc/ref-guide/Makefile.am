# can't change doxygen options via command line
# but we want to specify input and output directory
# to be compatible with out-of-source build
doxyinput = $(addprefix $(top_srcdir)/, src models BAT README.md)

define extradoxyconf
	PROJECT_NAME = "BAT"
	USE_MDFILE_AS_MAINPAGE = README.md
	GENERATE_LATEX = NO
	GENERATE_TAGFILE = bat-ref.tag
endef
export extradoxyconf

EXTRA_DIST = Doxyfile

include $(top_srcdir)/doc/Makefile-common.am

clean-local:
	-rm -rf html/ bat-ref.tag doxygen_sqlite3.db
