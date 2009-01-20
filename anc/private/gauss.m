function g = gauss(x, type, norm)
% g = gauss(x, type)
% Create a 1D gaussian.

if nargin < 2
    type = 'sigma';
end

if nargin < 3
    norm = 1;
end

if strfind(type, 'sigma')
    sig_x = x;
    lim_x = 3.*sig_x;
else
    lim_x = floor(x/2);
    sig_x = (lim_x - 2)/6;
end  

xx = -lim_x:lim_x;

if norm
    px = (sqrt(2*pi) * sig_x).^-1;
else
    px = 1;
end
g  = px .* exp(-0.5 .* (xx.^2 ./ sig_x.^2));
