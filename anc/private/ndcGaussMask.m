function [m, m_x, m_y, m_xy, m_xx, m_yy] = ndcGaussMask(s, type)
%function [m, m_x, m_y, m_xy, m_xx, m_yy] = ndcMask(s, type)
% 	Generate a mask for Normalised Convolution:

% Parse input arguments.
if nargin < 2
    type = 'dim';
end

if strfind(type, 'sigma') > 0
    lim = round(3*s);
elseif strfind(type, 'dim') > 0
    lim = s;
else
    error('Misunderstood type argument: %s', type);
end

if ~rem(lim, 2)
	lim = lim + 1;
end

[x, y] = meshgrid(-lim:lim, -lim:lim);

% Note: using gauss2 instead of knuttson.
m = gauss2(2*lim+1, 2*lim+1, 'dim');

m_x = m.*x;
m_y = m.*y;
m_xy = m.*(x.*y);
m_xx = m.*(x.^2);
m_yy = m.*(y.^2);
