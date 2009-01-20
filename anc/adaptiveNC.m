function [op, time, inform] = adaptiveNC(si, cm, varargin) 
% Interpolation using adaptive GST based normalised convolution (ANC)
% [op, time] = adaptiveNC(si, cm, varargin)
% 
% Required arguments:
%  si              = sampled image.
%  cm              = confidence map.
% Optional arguments:
%   (if left unspecified, defaults will be used):
%  'pad_amt'       = Amount to pad data.
%  'max_s'         = Maximum filter standard deviation.
%  'min_s'         = Minimum filter standard deviation.
%  'max_m'         = Maximum DoNC/NDC filter dimension.
%  'c'             = Standard deviation multiplier.
%  'alpha'         = Standard deviation exponent.
%  'anisoFiltSize' = Anisotropy filter standard deviation.
%  'gradFiltSize'  = Gradient filter standard deviation.        
%  'lamFiltSize'   = Eigenvalue filter standard deviation.        
%  'sigFiltSize'   = Standard deviation filter standard deviation(!). 
%  'phiFiltSize'   = Orientation filter standard deviation.  
%  'gradType'      = Type of filter to use for estimating gradients:
%                    'donc' or 'ndc'
%  'gradFiltType'  = DONC filter type:
%                    can be 'gauss', 'gauss_sig' or 'ndc'.
%  'verbose'       = Increase verbosity.
%  'type'          = What type of adaptive filtering to perform:
%                   'shape + scale' or 'scale'.
% Output Arguments:
%  op              = Output data.
%  time            = Time taken to process data.
%
% See Also:
%   griddata, kriging, rbf, natural_neighbour
%
% Matt Foster <ee1mpf@bath.ac.uk>

% Linted: 20060120 / 20090120

tic;

% Defaults:
opts = struct(	'pad_amt'				, 0,								...
								'max_s'					, 50,								...
								'min_s'					, 1.7,							...
								'max_m'					, 199,							...
								'c'							, 1,								...
								'alpha'					, 1.1,							...
								'anisoFiltSize' , 41,								...
								'gradFiltSize'	, 201,							...
								'lamFiltSize'		, 0,								... 
								'sigFiltSize'		, 0,								... 
								'phiFiltSize'		, 0,								...
								'gradFiltType'	, 'gauss',					...
								'verbose'				, 0,								...
								'gradType'			, 'donc',						... 
								'type'					, 'shape + scale',	...
								'padding'				, 'symmetric'       );

% Note: using struct here makes this function incompatible with octave
error(nargchk(2, 30, nargin, 'struct'));

% Parse the input arguments:
opts = args(varargin, opts, mfilename);

if sum(cm(:)) < 5
  warning('ANC:NotEnoughInputData', 'Input is very sparse');
  op	 = zeros(size(si));
  time = toc;
  return;
end

% Check for overridden pad_amt, and calculate if it hasn't been changed.
if ~opts.pad_amt
  opts.pad_amt = 7.* opts.max_s;
elseif opts.pad_amt < 6 * opts.max_s;
  warning('ANC:InsufficientPadding', 'Padding is probably too small.');
end

if opts.verbose
  fprintf(1, 'Padding: %d\n', opts.pad_amt); 
end

% Image size
[rows, cols] = size(si);

% Padding:
si = padarray(si, [opts.pad_amt, opts.pad_amt], opts.padding);
cm = padarray(cm, [opts.pad_amt, opts.pad_amt], opts.padding);

% Calculate bwdist:
sig_a = bwdist(cm); 
sig_a(sig_a == 0) = 1;

if strfind(opts.type, 'shape + scale')
  % Find filter size for donc.
  m_dist = ceil(max(max(sig_a)) + 2);
  if rem(m_dist, 2) == 0
    m_dist = m_dist + 5;
  end

  m_dist(m_dist > opts.max_m) = opts.max_m;

  if strfind(opts.gradType, 'donc')
    if opts.verbose; 
      fprintf(1, 'Running DoNC Edge detection. ');
      fprintf(1, 'Filter size: %d\n', m_dist); 
    end
    [gx, gy] = donc(si, cm, m_dist, opts.gradFiltType); 
  elseif strfind(opts.gradType, 'ndc')
    if opts.verbose 
      fprintf(1, 'Running NDC Edge detection. ');
      fprintf(1, 'Filter size: %d\n', m_dist); 
    end
    [gx, gy] = ndc(si, cm, m_dist, opts.gradFiltType);
  else
    error('Misunderstood edge detector type.');
  end

  if opts.verbose; 
    fprintf(1, 'Gradient Filter Size: %d.\n', opts.gradFiltSize); 
    fprintf(1, 'Generating gradient products.\n'); 
  end

  [gxy, gx2, gy2] = fastGenerateGradients(gx, gy, opts.gradFiltSize);

  if opts.verbose
    fprintf(1, 'Computing eigenvalues.\n'); 
  end

  [lam_1, lam_2] = gstEig(gxy, gx2, gy2);
  % Don't want to use the mex function? uncomment this:
  %[lam_1, lam_2] = gstEigenvalues(gxy, gx2, gy2);

  if opts.lamFiltSize
    if opts.verbose
      fprintf(1, 'Filtering Eigenvalues\n');
    end
    lf = gauss(opts.lamFiltSize, 'dim');
    lam_1 = sepConv2(lam_1, lf);
    lam_2 = sepConv2(lam_2, lf);
  end 

  if opts.verbose
    fprintf(1, 'Computing anisotropy.\n');
  end
  aniso = gstAnisotropy(lam_1, lam_2);

  if opts.anisoFiltSize
    if opts.verbose
      fprintf(1, 'Filtering anisotropy.\n'); 
    end
    af = gauss(opts.anisoFiltSize, 'dim');
    aniso = sepConv2(aniso, af);
  end

  if opts.verbose; 
    fprintf(1, 'Computing stdevs. c = %3.2f, alpha = %3.2f\n', ... 
    opts.c, opts.alpha);
  end

  [sig_u, sig_v] = gstSigmas(aniso, sig_a, opts.max_s, ...
  opts.min_s, opts.c, opts.alpha);

  clear aniso;

  if opts.sigFiltSize
    if opts.verbose
      fprintf(1, 'Filtering Sigma values.\n'); 
    end
    sf = gauss(opts.anisoFiltSize, 'dim');
    sig_u = sepConv2(sig_u, sf);
    sig_v = sepConv2(sig_v, sf);
  end

  if opts.verbose
    fprintf(1, 'Computing orientations.\n'); 
  end

  phi_2 = atan(gxy ./ (lam_2 - gy2));
  phi_2(isnan(phi_2)) = 0;

  if opts.phiFiltSize
    if opts.verbose 
      fprintf(1, 'Filtering Phi values.\n'); 
    end
    pf = gauss(opts.phiFiltSize, 'dim');
    phi_2 = sepConv2(phi_2, pf);
  end

  sig_u = sig_u(opts.pad_amt:opts.pad_amt+rows-1, ...
  opts.pad_amt:opts.pad_amt+cols-1);
  sig_v = sig_v(opts.pad_amt:opts.pad_amt+rows-1, ...
  opts.pad_amt:opts.pad_amt+cols-1);
  phi_2 = phi_2(opts.pad_amt:opts.pad_amt+rows-1, ...
  opts.pad_amt:opts.pad_amt+cols-1);
end

% Remove Padding:
sig_a = sig_a(opts.pad_amt:opts.pad_amt+rows-1, ...
opts.pad_amt:opts.pad_amt+cols-1);

% Construct filter arguments:
if strfind(opts.type, 'shape + scale')
	sig_x				= sig_v;
	sig_y				= sig_u;
	phi					= phi_2;
else
	sig_x				= sig_a;
	sig_y				= sig_x;
	phi					= zeros(size(sig_x)); % Cheat.
end

% Preallocate memory for outputs:
op_1 = zeros(rows, cols);
op_2 = zeros(rows, cols);

op = nan;

% check for nans
iterations = 1;
mult = 1:0.1:5;

while any(isnan(op(:))) && iterations <= length(mult)
	% Non-looping mex function.
	[op_1, op_2] = anisoDoubleFilter(double(si), double(cm), ...
	sig_y*mult(iterations), sig_x*mult(iterations), phi, opts.pad_amt);
	% Normalise the convolutions.
	op = op_1./op_2;
	iterations = iterations + 1;
end

inform.iterations = iterations;
inform.mult = mult(iterations-1); 

%op(isnan(op)) = op_2(isnan(op));	 

if any(isnan(op(:)))
	warning('ANC:GappyOutput', 'Output probably contains gaps')
end

time = toc;

end % of ANC

%----------------------------------------------------
% Utility subfunctions
%----------------------------------------------------
function [sig_u, sig_v] = gstSigmas(aniso, sig_a, max_s, min_s, c, alpha)
% function [sig_u, sig_v] = gstSigmas(aniso, sig_a, min_s, c, alpha)
%		Compute gst standard deviations.

if nargin < 5
  c = 0.5;
end
if nargin < 6
  alpha = 0.5;
end

sig_u = c.*(1+aniso).^alpha.*sig_a;
sig_v = c.*(1-aniso).^alpha.*sig_a;

sig_u(isnan(sig_u)) = 0;
sig_v(isnan(sig_v)) = 0;
sig_u(sig_u > max_s) = max_s;
sig_v(sig_v > max_s) = max_s;
sig_u(sig_u < min_s) = min_s;
sig_v(sig_v < min_s) = min_s;

end % gstSigmas

function [aniso, energy] = gstAnisotropy(lam_1, lam_2)
  % function [aniso, energy] = gstAnisotropy(lam_1, lam_2)
  % Compute anisotropy and energy from given eigenvaulues.

  energy = lam_1 + lam_2;
  diff = lam_1 - lam_2;
  aniso = diff ./ energy;

  aniso(energy == 0) = 0; % necessary?
  aniso(diff == 0) = 0;

end % gstAnisotropy

% Only used when mex is commented. 
function [lam_1, lam_2] = gstEigenvalues(gxy, gx2, gy2)
  % function [lam_1, lam_2] = gstEigenValues(gxy, gx2, gy2)
  % Create gst eigenvalues, and sort them into lam_1 and lam_2.

  [rr, cc] = size(gxy);

  % Preallocation for speed:
  lam_1 = zeros(size(gxy));
  lam_2 = zeros(size(gxy));

  % Loop, calculating the eigenvalues.
  for i = 1:rr
    for j = 1:cc
      eigenValues = eig([gx2(i,j), gxy(i,j); gxy(i,j), gy2(i,j)]);
      lam_1(i,j) = max(eigenValues);
      lam_2(i,j) = min(eigenValues);
    end
  end

end % gstEigenvlaues

function opts = args(varargin, opts, scriptname)
  % opts = args(varargin, opts, scriptname)
  % Parse input arguments into a struct.

  if mod(length(varargin),2)
    error('Arguments should be in pairs.');
  end

  % Argument handler.
  for i = 1:2:length(varargin)
    name = varargin{i};
    value = varargin{i+1};

    ccmd = sprintf('cls = class(opts.%s);', name);
    eval(ccmd);
    if ~strcmpi(class(value), cls)
      warning('ANC:Arguments:IncorrectClass', ...
        'Wrong class in option %s. Expected %s', ...
        name, cls);
    end

    % Nice way of accessing struct elements.
    opts.(name) = value;
  end

  if isfield(opts, 'verbose')
    if opts.verbose > 0
      fprintf(1, '\n%s options:\n\n', scriptname);
      disp(opts)
    end
  end

end % args

