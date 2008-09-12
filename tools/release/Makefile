
DIRBASE   = bat
VERSION   = 0.1

RELDIR0   = $(DIRBASE)-$(VERSION)

SOURCEDIR = /home/dkollar/work/BayesianAnalysisToolkit/BAT

DESTDIR   = .

CVSCLEANER = $(shell pwd)/removeCVS.sh

#SOURCEFILES := $(wildcard $(SOURCEDIR)/src/*.cxx)
#INCFILES    := $(wildcard $(SOURCEDIR)/inc/*.h)
#OTHERINCS   := $(wildcard $(SOURCEDIR)/inc/*.c)

REL_DIRS_ST = src include
REL_DIRS_EMPTY = obj

DIRS_ALL  = $(REL_DIRS_ST) $(REL_DIRS_EMPTY) examples doc

FILES_ALL = Makefile

################
# targets below

all: tar


credir:
	mkdir $(RELDIR0)

copyall: credir files $(DIRS_ALL) version_header
	

files:
	cp $(SOURCEDIR)/Makefile $(RELDIR0)

$(REL_DIRS_ST):
	cp -r $(SOURCEDIR)/$@ $(RELDIR0)
	cd $(RELDIR0)/$@; $(CVSCLEANER)

examples:
	cp -r $(SOURCEDIR)/$@ $(RELDIR0)
	cd $(RELDIR0)/$@; $(CVSCLEANER) example05; $(CVSCLEANER) example06; $(CVSCLEANER)

$(REL_DIRS_EMPTY):
	mkdir $(RELDIR0)/$@

tools:
	cp -r $(SOURCEDIR)/$@ $(RELDIR0)
	cd $(RELDIR0)/$@; $(CVSCLEANER);  $(CVSCLEANER) analyzeMCMC
	rm $(RELDIR0)/$@/setup.sh

doc: introduction.tex
	mkdir $(RELDIR0)/$@
	sed -e "s/__CVS_VERSION__/$(VERSION)/" $(SOURCEDIR)/$@/README > $(RELDIR0)/README
	cp $(SOURCEDIR)/$@/INSTALL $(RELDIR0)
	cp $(SOURCEDIR)/$@/COPYING $(RELDIR0)/$@
	cp $(SOURCEDIR)/$@/LICENSE $(RELDIR0)/$@
	cp $(SOURCEDIR)/$@/releasenotes.txt $(RELDIR0)/$@
	cp $(SOURCEDIR)/$@/changelog $(RELDIR0)/$@
	mv introduction.ps $(RELDIR0)/$@

introduction.tex:
	mkdir introduction-tmp
	cp -r $(SOURCEDIR)/doc/introduction introduction-tmp
	cd introduction-tmp/introduction; latex introduction; latex introduction; dvips introduction
	mv introduction-tmp/introduction/introduction.ps .
	rm -rf introduction-tmp

version_header: include
	sed -e "s/cvs/$(VERSION)/"  $(RELDIR0)/include/BCReleaseVersion.h > BCReleaseVersion.h.ch
	mv BCReleaseVersion.h.ch $(RELDIR0)/include/BCReleaseVersion.h

tar: copyall
	tar -czf $(RELDIR0).tgz $(RELDIR0)



