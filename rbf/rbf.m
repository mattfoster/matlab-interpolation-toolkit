function [varargout] = rbf(xi, yi, zi, xx, yy, basis_func, lambda, log_fudge)
% RBF Perform radial basis function interpolation.
% zz = rbf(xi, yi, zi, x, y, basis_function)
%
% Interpolate the scattered values xi, yi, zi at xx, yy (which should be plaid)
%
% Functionality is similar to griddata.
%
% Available basis functions:
% euclidean                             Default. Also fits first order
%                                       polynomial surface.
% thin_plate_spline / biharmonics       These are the same.
%                                       Also fits polynomial surface.
% gaussian                              Constant not implemented.
% multiquadratic                        Constant not implemented.
% triharmonic
%
% lambda:
% regularisation parameter see [2]. With regularisation, the interpolation
% condition is relaxed, and the fitted surface will be smoothed. Setting this
% to zero implies interpolation, and high values reduce to a least-squares
% affine model.
%
% log_fudge: 
% what to do when when radius = 0. Generally you'll want this to be 0 too.
%
% Example:
% [yy, xx]     = meshgrid(1:100, 1:100);
% % Make some correlated data.
% input        = corr_data(100, 15);
% aa           = rand(100) > 0.95;
% sampled      = input.*aa;
% [xi, yi, zi] = find(sampled);
%
% [zz, time]   = rbf(xi, yi, zi, xx, yy);
% [zz_tps, time] = rbf(xi, yi, zi, xx, yy, 'thin_plate_spline')
%
% rmse_zz     = sqrt(mean((input(:) - zz(:)).^2))
% rmse_zz_tps = sqrt(mean((input(:) - zz_tps(:)).^2))
%
%
% Reference: 
% [1] Surface interpolation with radial basis functions, 
%       Carr et al, 
%       IEEE Trans. on Medical Imaging.
%       Vol 16. No. 1, Feb. 1997
%
% [2] Shape matching and object recognition using contexts.
%       Belongie et al,
%       IEEE Trans. Pattern Analysis and Machine Intelligence
%       Vol. 24, No. 24, April 2002
%
% Matt Foster <ee1mpf@bath.ac.uk>

% Start timer
tic;

% Check input and output argument numbers
error(nargchk(5, 8, nargin, 'struct'));
error(nargoutchk(1, 4, nargout, 'struct'))

% Change this to control the behaviour of log(0). 
if nargin < 8
  log_fudge = 0;
end

if nargin >= 6
  basis_func = str2func(basis_func);
else
  basis_func = @euclidean;
end

[x1,x2]      = meshgrid(yi);
[y1,y2]      = meshgrid(xi);

[basis, poly_order] = basis_func(x1, x2, y1, y2, log_fudge);

if nargin >= 7
  scale = euclidean(x1, x2, y1, y2);
  scale = mean(scale(:));
  lambda = scale.^2 .* lambda;
else
  lambda = 0;
end

% Perform regularisation (or not)
basis = basis + lambda * eye(size(basis));

% So far I've only seen first order polynomials, or none
if poly_order == 0
  mat = basis;
  ff  = zi;
elseif poly_order == 1 
  q   = [ones(length(xi), 1), xi, yi];
  mat = [basis, q; q', zeros(3)];
  ff  = [zi; zeros(3,1)];
end

% Things below could get quite memory intensive.
clear basis;

% Get the coefficients.
lam_c = mat \ ff;

% Extract the poly coefficients.
lam   = lam_c(1:end-poly_order*2-1*poly_order);
c     = lam_c(end-poly_order*2:end);
clear lam_c;
clear mat;

% Now reconstruct everything:
if poly_order > 0
  poly = c(1) +  xx.*c(2) + c(3).*yy;
else
  poly = zeros(size(xx));
end

% Create matric for output
zz = zeros(size(xx));

%----------------------------------------------------
% for ii = 1:length(xx(:))
%   % loop over inputs
%   for kk = 1:length(zi)
%     norm = basis_func(xx(ii), xi(kk), yy(ii), yi(kk), log_fudge);
%     zz(ii) = zz(ii) + lam(kk).*norm;
%   end
% end
%----------------------------------------------------

% This is a vectorised version of the above:
for ii = 1:length(xx(:))

  nm = basis_func( ...
    repmat(xx(ii), [size(xi(:), 1), 1]), ...
    xi(:), ... 
    repmat(yy(ii), [size(yi(:), 1), 1]), ...
    yi(:), ...
    log_fudge);

  zz(ii) = sum(lam .* nm);
end

% Add the polynomial
zz = zz + poly;

% stop timer
time = toc;

% 4 output args -> x, y, z, time
% 3 output args -> x, y, z
% 2 output args -> z, time
% 1 output arg  -> z
if nargout == 4
  varargout{1} = xx;
  varargout{2} = yy;
  varargout{3} = zz;
  varargout{4} = time;
elseif nargout == 3
  varargout{1} = xx;
  varargout{2} = yy;
  varargout{3} = zz;
elseif nargout == 2
  varargout{1} = zz;
  varargout{2} = time;
elseif nargout == 1 
  varargout{1} = zz;
end

end % rbf 

% Basis functions
% mlint says these are never used.. it's wrong!


% Poly is an integer detailing the order of polynomial to use.

% Euclidean distance -- Linear basis function
function [bf, poly] = euclidean(x1, x2, y1, y2, varargin) %#ok
  bf = sqrt((x1 - x2).^2 + (y1 - y2).^2);
  poly = 1;
end

% Thin plate spline 
function [bf, poly] = thin_plate_spline(x1, x2, y1, y2, varargin) %#ok
  rr = euclidean(x1, x2, y1, y2);
  % Not 100% sure what base to use for the log.
  bf = rr.^2 .* log(rr);
  bf(rr == 0) = varargin{1};
  poly = 1;
end

% Gaussian
function [bf, poly] = gaussian(x1, x2, y1, y2, varargin) %#ok
  rr = euclidean(x1, x2, y1, y2);
  if nargin < 5
    aa = 1;
  else
    aa = varargin{1};
  end
  bf = exp(aa.*rr.^2);
  % No polynomial here.
  poly = 0;
end

% Multiquadratic
function [bf, poly] = multiquadratic(x1, x2, y1, y2, varargin) %#ok
  rr = euclidean(x1, x2, y1, y2);
  if nargin < 5
    cc = 1;
  else
    cc = varargin{1};
  end
  bf = sqrt(rr.^2 + cc.^2);
  poly = 0;
end

% Triharmonic Spline -- C2 Continuity
function [bf, poly] = triharmonic_spline(x1, x2, y1, y2, varargin) %#ok
  rr = euclidean(x1, x2, y1, y2);
  bf = rr.^4 .* log(rr);
  bf(rr == 0) = varargin{1};
  poly = 0;
end
