function op = sepConv2(in, filt)
% op = sepConv2(in, fs, type)
% Convolve a matrix with a 2D gaussian using a separable transform.
%
op = conv2(in,  filt, 'same');
op = conv2(op,  filt', 'same');
