# can't change doxygen options via command line
# but we want to specify input and output directory
# to be compatible with out-of-source build
doxyinput = $(srcdir)/front.md
doxyinput_raw = src models BAT
doxyinput += $(doxyinput_raw:%=$(top_srcdir)/%)

define extradoxyconf
	PROJECT_NAME = "BAT"
	USE_MDFILE_AS_MAINPAGE = front.md
	IMAGE_PATH = $(srcdir)/..
	GENERATE_LATEX = NO
	GENERATE_TAGFILE = bat-ref.tag
endef
export extradoxyconf

# do nothing by default so `make` in top-level dir doesn't required all the
# tools used for the docs
all: ;

include $(top_srcdir)/doc/Makefile-common.am

# don't rerun doxygen if file exists. `all` always does that
bat-ref.tag:
	$(MAKE) doxy

clean-local:
	-rm -rf html/ bat-ref.tag doxygen_sqlite3.db
