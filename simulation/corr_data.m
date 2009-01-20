function c = corr_data(dim, disk)
% Create a random signal, and impose correlation.
% c = corr_data(dim, disk)  
% 
% Required arguments:
%   'dim' = size of the field.
%   'disk' = size of the disc used to filter the field.
%
% Literature: 
%   THE VARIOGRAM AND ITS ESTIMATION, Omre, 1984 
%
% See also:
%   corr_data_uni
%
% Matt Foster <ee1mpf@bath.ac.uk>

if nargin < 2
    error('Expected 2 arguments')
end

if any(size(disk) ~= 1)
    error('Disc size should be a single integer')
end

a = randn(dim); 
f = fspecial('disk', disk);
c = filter2(f, a);

c = c + abs(min(c(:)));
c = c./max(c(:));
