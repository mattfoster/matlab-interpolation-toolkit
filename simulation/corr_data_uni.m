function c = corr_data_uni(dim, disk, pad, nh)
% Create a random signal, and impose correlation (lognormal version).
% c = corr_data_uni(dim, disk, pad, nh)  
% 
% Required arguments:
%   'dim' = size of the field.
%   'disk' = size of the disc used to filter the field.
% Optional arguments:
%   'pad' = size of padding (default 50)
%   'nh'  = size of neighbourhood (default 2) used for multinormality reduction.
%
% Literature: 
%   THE VARIOGRAM AND ITS ESTIMATION, Omre, 1984 
%
% See also:
%   corr_data
%
% Matt Foster <ee1mpf@bath.ac.uk>

error(nargchk(2, 4, nargin));

if nargin < 3
    pad = 50;
end

if nargin < 4
    nh = 2; % 2*nh+1 square neighbourhood
end

if any(size(disk) ~= 1)
    error('Disc size should be a single integer')
end


% Get some multinormal data
d = corr_data(dim+(2*pad), disk);

% remove multinormality
% Look at local neighbourhood 5x5 mask choose randomly from
% top 10

for ii = pad:size(d,1)-pad-1
    for jj = pad:size(d,nh)-pad-1
        nhood = d(ii-nh:ii+nh, jj-nh:jj+nh);
        sv = sort(nhood(:),'descend');
        pos = ceil(10*rand(1));
        new(ii-pad+1,jj-pad+1) = sv(pos);
    end
end

% Make the histogram look more like omre's
c = abs(reallog(new));
