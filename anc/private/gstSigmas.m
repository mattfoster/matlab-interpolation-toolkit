function [sig_u, sig_v] = gstSigmas(aniso, sig_a, max_s, min_s, c, alpha)
% function [sig_u, sig_v] = gstSigmas(aniso, sig_a, min_s, c, alpha)
% 	Compute gst standard deviations.

if nargin < 5
	c = 0.5;
end
if nargin < 6
	alpha = 0.5;
end

sig_u = c.*(1+aniso).^alpha.*sig_a;
sig_v = c.*(1-aniso).^alpha.*sig_a;

%----------------------------------------------------
% sig_u = sig_a .* alpha ./ (alpha + aniso);
% sig_v = sig_a .* (alpha + aniso) ./ aniso;
%----------------------------------------------------

sig_u(isnan(sig_u)) = 0;
sig_v(isnan(sig_v)) = 0;
sig_u(sig_u > max_s) = max_s;
sig_v(sig_v > max_s) = max_s;
sig_u(sig_u < min_s) = min_s;
sig_v(sig_v < min_s) = min_s;
