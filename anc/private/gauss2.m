function g = gauss2(x, y, type, norm)
% g = gauss2(x, y, type)
% Create a 2D gaussian by convolution.

if nargin < 3
    type = 'sigma';
end

if strfind(type, 'sigma')
    sig_x = x;
    sig_y = y;
    lim_x = 3.*sig_x;
    lim_y = 3.*sig_y;
else
    lim_x = floor(x/2);
    lim_y = floor(y/2);
    sig_x = (lim_x - 2)/6;
    sig_y = (lim_y - 2)/6;
end  

xx = -lim_x:lim_x;
yy = -lim_y:lim_y;

px = (sqrt(2*pi) * sig_x).^-1;
fx = px .* exp(-0.5 .* (xx.^2 ./ sig_x.^2));

py = (sqrt(2*pi) * sig_y).^-1;
fy = py .* exp(-0.5 .* (yy.^2 ./ sig_y.^2));

g = conv2(fx, fy');
