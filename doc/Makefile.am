SUBDIRS = ref-guide manual

EXTRA_DIST = \
	COPYING CREDITS LICENSE releasenotes.md \
	Doxyfile header.html footer.html \
	customdoxygen.css doxy-boot.js newcommands.js

.PHONY: all
all: manual

.PHONY: ref-guide
ref-guide:
	make -C ref-guide/ doxy

# create html and pdf via latex
# first need the ref-guide for its tag file
.PHONY: manual
manual: ref-guide
	make -C manual/ pdf
