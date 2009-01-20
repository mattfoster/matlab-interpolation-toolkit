##
## Mex Options:
##

# Will probably need customising for your system

CC         = gcc
MEX       = /Applications/MATLAB_R2008b.app/bin/mex
DEBUG      = -g -O -v
MEXSUFFIX  = $(/Applications/MATLAB_R2008b.app/bin/mexext)

##
## Linker/Loader flags
##

LDFLAGS    =
INCLUDES   = -lm
CFLAGS     = -g $(INCLUDES)
PROGRAM	   = anisoDoubleFilter
FUNCTION   = $(PROGRAM)$(MEXSUFFIX)
MAKEFILE   = Makefile
HDRS	   = defs.h 
SRCS	   = anisoDouble.c common.c
MAINSRC	   = $(PROGRAM).c

# Script to generate headers.
#HEADER     = include.sh # requires cproto.
LIBS       = $(LOCAL_LIBS) $(SYS_LIBS)

##
## Targets
##

#all: include $(FUNCTION) #standalone
all: $(FUNCTION) #standalone

$(FUNCTION): $(MAKEFILE)

sources : $(SRCS) $(HDRS)

include: $(SRCS) 
	@echo 'Generating headers for $+'
	$(HEADER) $+

test : 	$(TESTSCRIPT) $(FUNCTION)
	cat $(TESTSCRIPT)
	$(MATLAB)/bin/matlab < $(TESTSCRIPT)

$(FUNCTION): $(MAINSRC) $(SRCS) $(HDRS) Makefile
	$(MEX) $(DEBUG) $(CFLAGS) $(LDFLAGS) $(SRCS) $(MAINSRC) \
	$(LIBS) -o $(FUNCTION)

standalone: main.c $(SRCS) 
	$(CC) -o anisoGauss -Wall -ansi $(CFLAGS) $(LDFLAGS) $+

delmex:
	rm -f $(FUNCTION)

clean : delmex 
	rm -f *~ core

tags : sources
	etags -t $(SRCS) $(HDRS)

.PHONY:    all sources test delmex full
.PHONY:    clean tags
