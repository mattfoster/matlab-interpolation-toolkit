function zi = natural_neighbour(xx, yy, zz, xi, yi)
% Perform interpolation using the natural neighbour method.
% zi = natural_neighbour(xx, yy, zz, xi, yi)
%
%	Perform natural neighbour interpolation using code from 
% http://www.marine.csiro.au/~sak007/
% 
% Arguments:
%   xx, yy, zz = scattered input data 
%   xi, yi     = output positions
%
% See also:
%   griddata, adaptiveNC, rbf, kriging
%
% Mex gateway by Matt P. Foster <ee1mpf@bath.ac.uk>

error(nargchk(5, 5, nargin));
zi = mex_nn(xx, yy, zz, xi, yi);

end % function