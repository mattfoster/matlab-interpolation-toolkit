# Makefile for Matlab interpolation kit

# This may need setting manually
# EXT = /usr/local/matlab2006b/bin/mexext
# EXT = /Applications/MATLAB_R2008b.app/bin/mexext

all: anc natural_neighbour toolkit contents

anc: 
	make -C anc/mex/gst_eigenvalue all
	make -C anc/mex/anisotropic_filter all

natural_neighbour:
	$(MAKE) -C natural_neighbour

toolkit:
	mkdir -p toolkit/private
	find . -path './toolkit' -prune -o -name '*.mex*' -exec cp {} toolkit/private \;
	find . -path './toolkit' -prune -o -maxdepth 2 -name '*.m' -exec cp {} toolkit \;
	cp anc/private/*.m toolkit/private/

contents: toolkit
	python contents.py
	
clean:
	$(MAKE) -C anc/mex/gst_eigenvalue clean
	$(MAKE) -C anc/mex/anisotropic_filter clean
	$(MAKE) -C natural_neighbour clean
	rm -r toolkit

.PHONY: clean anc natural_neighbour toolkit
