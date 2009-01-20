function [gxy, gx2, gy2] = fastGenerateGradients(gx, gy, filtSize);
% [gxy, gx2, gy2] = fastGenerateGradients(gx, gy, filtSize)
% Find gradiants using seperable filtering.
%
% Matt Foster <ee1mpf@bath.ac.uk>

% $LastChangedDate$

gxy = gx.*gy;
gx2 = gx.^2;
gy2 = gy.^2;

gf = gauss(filtSize, 'dim'); % 1D gaussian

gxy = sepConv2(gxy, gf);
gx2 = sepConv2(gx2, gf);
gy2 = sepConv2(gy2, gf);
