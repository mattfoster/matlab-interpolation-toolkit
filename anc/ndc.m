function [gx, gy] = ndc(si, cm, r, type, al, be)
% Calculate the derivative of a sampled image using NDC.
% [xx, yy] = ndc(si, cm, al, be, r)
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
%   ADAPTIVENC, DONC
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
    al = 3;
    be = 0.5;
  end
  % Create a Knuttson mask
  [m, m_x, m_y, m_xy, m_xx, m_yy] = ndcMask(al, be, r);
elseif strfind(type, 'gauss') > 0
  if strfind(type, 'gauss_sig') > 0
    g_type = 'sigma';
  else
    g_type = 'dim';
  end
  % Create a guassian mask
  [m, m_x, m_y, m_xy, m_xx, m_yy] = ndcGaussMask(r, g_type);
  
else
  error('Misunderstood type argument: %s', type);
end

cm 	  = double(cm); 

nc    = conv2(cm, m, 'same');
cc    = conv2(si, m, 'same');

c_x   = conv2(si, m_x, 'same');
c_y   = conv2(si, m_y, 'same');

nc_x  = conv2(cm, m_x, 'same');
nc_y  = conv2(cm, m_y, 'same');

d_x	  = nc .* c_x - nc_x .* cc;
d_y	  = nc .* c_y - nc_y .* cc;

n_xx  = nc .* conv2(cm, m_xx, 'same') - nc_x.^2;
n_xy  = nc .* conv2(cm, m_xy, 'same') - nc_x .* nc_y;
n_yx  = n_xy; 
n_yy  = nc .* conv2(cm, m_yy, 'same') - nc_y.^2;

[rows, cols] = size(si);
for i = 1:rows
	for j = 1:cols
		D = [d_x(i, j); d_y(i, j)];
		N = [n_xx(i, j), n_xy(i, j); n_yx(i, j), n_yy(i, j)];
		G = inv(N)*D;
		gx(i, j) = G(1, 1);
		gy(i, j) = G(2, 1);
	end
end


