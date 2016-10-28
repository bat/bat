SUBDIRS = ref-guide manual

EXTRA_DIST = \
	COPYING CREDITS LICENSE releasenotes.md \
	customdoxygen.css  newcommands.js

.PHONY: ref-guide
ref-guide:
	make -C ref-guide/ doxy

# create html and pdf via latex
# first need the ref-guide for its tag file
.PHONY: manual
manual: ref-guide
	make -C manual/ pdf
