function [gx, gy] = donc(si, cm, r, type, al, be)
% Calculate the derivative of a sampled image using DoNC.
% [gx, gy] = donc(si, cm, r, type, al, be)
% DoNC has a lower complexity than plain NDC, and a slightly lower output 
% quality.
% 
% Required arguments:
%   si          = sampled image.
%   cm          = confidence map.
% Optional arguments:
%   r           = mask radius (default 11).  
%   'type'      = set the filter type:
%     'gauss'   = gaussian.
%     'kuttson' = Knuttson style raised cosine.
%   'al'        = Knuttson filter alpha param.
%   'be'        = Knuttson filter beta param.
%
% See also:
%   ADAPTIVENC, NDC
%   
% Matt Foster <ee1mpf@bath.ac.uk>

error(nargchk(2, 6, nargin));

if nargin < 3
  r = 11;
end

if nargin < 4
  type = 'gauss';
end

if strfind(type, 'knuttson') > 0
  if nargin < 5
    al = 0;
    be = 2;
  end
  % Create a Knuttson mask
  [m, m_x, m_y] = ndcMask(al, be, r);
elseif strfind(type, 'gauss') > 0
  if strfind(type, 'gauss_sig') > 0
    g_type = 'sigma';
  else
    g_type = 'dim';
  end
  % Create a guassian mask
  [m, m_x, m_y] = ndcGaussMask(r, g_type);
else
  error('Misunderstood type argument: %s', type);
end

cm    = double(cm); 

nc    = conv2(cm, m, 'same');
cc    = conv2(si, m, 'same');

c_x   = conv2(si, m_x, 'same');
c_y   = conv2(si, m_y, 'same');

nc_x  = conv2(cm, m_x, 'same');
nc_y  = conv2(cm, m_y, 'same');

% Not actually used:
%----------------------------------------------------
% d_x   = nc .* c_x - nc_x .* cc;
% d_y   = nc .* c_y - nc_y .* cc;
%----------------------------------------------------

gx   = (c_x .* nc - nc_x .* cc) ./ (nc.^2);
gy   = (c_y .* nc - nc_y .* cc) ./ (nc.^2);

gx(isnan(gx)) = 0;
gy(isnan(gy)) = 0;
