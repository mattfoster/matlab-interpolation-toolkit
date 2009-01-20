#!/usr/bin/env python

from __future__ import with_statement

import glob
import os.path
import sys

header = """\
% Matlab Scattered Data Interpolation Toolbox
%
% Author:  Matt Foster <ee1mpf@bath.ac.uk>, CSAOS, University of Bath
% License: BSD (Except: natural neighbour -- see readme for license)
%
% Interpolation Functions:"""

footer = """\
%
% Demos:
%
% To generate correlated test data, run:
%
% [yy, xx]     = meshgrid(1:100, 1:100);
% % Make some correlated data.
% input        = corr_data(100, 15);
% cm           = rand(100) > 0.95;
% sampled      = input.*cm;
% [xi, yi, zi] = find(sampled);
%
% figure
% imagesc(sampled)
%
% This will generate some sparse test data suitable for both adaptiveNC and
% griddata style methods. To interpolate this using ANC, run:
%
% anc_out = adaptiveNC(sampled, cm)
% figure
% imagesc(anc_out)
%
% To interpolate the data using griddata style functions (RBF, kriging, 
% natural neighbour, etc.), run:
%
% nat_out = natural_neighbour(xx, yy, zz, xx, yy);
% rbf_out = rbf(xx, yy, zz, xx, yy);
% kriging_out = kriging(xx, yy, zz, xx, yy);
% figure
% imagesc(nat_out) % etc.
%
% Finally, calculate RMSE values, using:
% anc_rmse     = sqrt(mean((input(:) - anc_out(:)).^2))
% nat_rmse     = sqrt(mean((input(:) - nat_out(:)).^2))
% % etc.
%
% Please note: Kriging and RBF will probably not work on large data sets due to 
% large matrix inversions. You mileage may vary.
% """

with open('toolkit/Contents.m', 'w') as output:
    sys.stdout = output # redirect stdout to output file
    print header
    for fname in glob.glob('toolkit/*.m'):
        # Skip this iteration if we're reading the output file
        if fname.endswith('Contents.m'):
            continue
        with open(fname,'r') as f:
            f.readline()        # read past the first line
            line = f.readline() # read the second
        base = os.path.basename(fname)
        name = os.path.splitext(base)[0]
        print "%% %-20s\t- %s" % (name, line.strip('%').rstrip())
    print footer
    
