MATLAB Multivariate Interpolation Toolbox
=========================================

This toolbox contains code for 2D multivariate interpolation in MATLAB.

Contents
========

* Adaptive Normalised Convolution (ANC)
* Radial Basis Function Interpolation (RBF)
* Kriging
* Natural Neighbour Interpolation
  * The included Natural Neighbours interpolation software was provided by
    Pavel Sakov and is available at URL:
    http://www.marine.csiro.au/~sakov/nn.tar.gz

Unless otherwise stated, all code is Copyright Matt Foster <ee1mpg@bath.ac.uk>

Installation
============

Some functions require MEX functions to built. To compile all of the necessary
Makefiles in UNIX, type `make`. Some of the Makefiles may require alterations
for your system, however, provided you have a complete build environment, and
the binaries `mex` and `mexext` are available, the process should be fairly
simple. You should then find a directory named `toolkit`. Copy this directory somewhere in your MATLAB work area, and finally, from within MATLAB, issue the command:

    addpath('interpolation_toolkit');

It is not necessary to use `genpath`, or include the `private` directory

Usage
=====

With the exception of ANC, which interpolates onto matrices, all of the
functions provided have prototypes of the following form:

    [xx, yy, zz] = function(xi, yi, zi, xx, yy, [optional args])

or more simply
    
    zz = function(xi, yi, zi, xx, yy, [optional args])

The help within the various functions gives more information on the input and
output arguments. For example, issuing the command:

    help rbf

Will print the built in help for the RBF interpolation command.

Running the command:

    help interpolation_toolkit

Where `interpolation_toolkit` is the directory name, will give basic help for the toolkit, including some usage instructions.

ANC
===

ANC operates using matrices, and requires arguments of the following form:

    out = adaptiveNC(si, cm, [optional args])

Full help is available within the m-file.

Bugs
====

Please report bugs to Matt Foster <ee1mpf@bath.ac.uk>

History
=======

Moved from local SVN to GitHub: 20090120
