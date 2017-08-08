function b = rcosdesign(beta, span, sps, varargin)
%RCOSDESIGN Raised cosine FIR filter design
%   B = rcosdesign(BETA, SPAN, SPS) returns square root raised cosine FIR
%   filter coefficients, B, with a rolloff factor of BETA. The filter is
%   truncated to SPAN symbols and each symbol is represented by SPS
%   samples. RCOSDESIGN designs a symmetric filter. Therefore, the filter
%   order, which is SPS*SPAN, must be even. The filter energy is one.
%
%   B = rcosdesign(BETA, SPAN, SPS, SHAPE) returns a normal raised cosine
%   FIR filter when you set SHAPE to 'normal'. When you set SHAPE to
%   'sqrt', the function returns a square root raised cosine filter.
%   
%   % Example 1:
%   %   Design a square root raised cosine filter with a rolloff factor of
%   %   0.25. Truncate the filter to 6 symbols and represent each symbol
%   %   with 4 samples.
%
%   h = rcosdesign(0.25, 6, 4);        % Raised cosine FIR filter design
%   fvtool(h, 'Analysis', 'impulse')   % Visualize the filter
%
%   % Example 2:
%   %   Interpolate using a square root raised cosine filter with a rolloff 
%   %   factor of 0.25, and a filter span of 6. Use an interpolation factor
%   %   of 4. Note that the upfirdn function flushes out all the input 
%   %   samples by appending zeros at the end of the input.
%
%   sps = 4;
%   h = rcosdesign(0.25, 6, sps);      % Raised cosine FIR filter design
%   d = randi([0 1], 100, 1);
%   x = upfirdn(d, h, sps);
%   plot(x)
%
%   See also gaussdesign.

%   Copyright 2012-2013 The MathWorks, Inc.

% Argument error checking
narginchk(3,4)

validateattributes(beta, {'double', 'single'}, ...
  {'real', 'nonnegative', '<=', 1, 'scalar'}, ...
  '', 'BETA', 1)

if beta == 0,
   beta = realmin;
end

validateattributes(span, {'double', 'single'}, ...
  {'real', 'positive', 'scalar', 'finite'}, ...
  '', 'SPAN', 2)

validateattributes(sps, {'double', 'single'}, ...
  {'real', 'positive', 'scalar', 'finite', 'integer'}, ...
  '', 'SPS', 3)

if mod(sps*span, 2) == 1
  error(message('signal:rcosdesign:OddFilterOrder', sps, span))
end

if nargin > 3
  shape = varargin{1};
else
  shape = 'sqrt';
end

shape = validatestring(shape, {'sqrt', 'normal'}, '', 'SHAPE', 4);


% Design the raised cosine filter

delay = span*sps/2;
t = (-delay:delay)/sps;

if strncmp(shape, 'normal', 1)
  % Design a normal raised cosine filter
  
  % Find non-zero denominator indices
  denom = (1-(2*beta*t).^2);
  idx1 = find(abs(denom) > sqrt(eps));
  
  % Calculate filter response for non-zero denominator indices
  b(idx1) = sinc(t(idx1)).*(cos(pi*beta*t(idx1))./denom(idx1))/sps;
  
  % fill in the zeros denominator indices
  idx2 = 1:length(t);
  idx2(idx1) = [];
  
  b(idx2) = beta * sin(pi/(2*beta)) / (2*sps);
  
else
  % Design a square root raised cosine filter
  
  % Find mid-point
  idx1 = find(t == 0);
  if ~isempty(idx1),
    b(idx1) = -1 ./ (pi.*sps) .* (pi.*(beta-1) - 4.*beta );
  end
  
  % Find non-zero denominator indices
  idx2 = find(abs(abs(4.*beta.*t) - 1.0) < sqrt(eps));
  if ~isempty(idx2),
    b(idx2) = 1 ./ (2.*pi.*sps) ...
      * (    pi.*(beta+1)  .* sin(pi.*(beta+1)./(4.*beta)) ...
      - 4.*beta     .* sin(pi.*(beta-1)./(4.*beta)) ...
      + pi.*(beta-1)  .* cos(pi.*(beta-1)./(4.*beta)) ...
      );
  end
  
  % fill in the zeros denominator indices
  ind = 1:length(t);
  ind([idx1 idx2]) = [];
  nind = t(ind);
  
  b(ind) = -4.*beta./sps .* ( cos((1+beta).*pi.*nind) + ...
    sin((1-beta).*pi.*nind) ./ (4.*beta.*nind) ) ...
    ./ (pi .* ((4.*beta.*nind).^2 - 1));
  
end

% Normalize filter energy
b = b / sqrt(sum(b.^2));
end
