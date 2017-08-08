function [DH,DW] = highpass(N, F, GF, W, delay)
%HIGHPASS Desired frequency response for highpass filters.
%  CFIRPM(N,F,@highpass, ...) designs a linear-phase highpass filter
%  response using CFIRPM.
%
%  CFIRPM(N,F,{@highpass, D}, ...) specifies group-delay offset D such
%  that the filter response will have a group delay of N/2 + D in units of
%  the sample interval, where N is the filter order. Negative values create
%  less delay, while positive values create more delay.  By default, D=0.
%
%  The symmetry option SYM defaults to 'even' if unspecified in the call to
%  CFIRPM, if no negative band edge frequencies are specified in F.
%
%  EXAMPLE: Design a 31-tap, complex highpass filter.
%    b = cfirpm(30,[-1 -.5 -.4 .7 .8 1],@highpass);
%    fvtool(b);  % View filter response.
%
%  EXAMPLE: Reduced group delay filter response.
%    b = cfirpm(30,[0 .6 .7 1],{@highpass,-1});
%    fvtool(b);
%
%  See also CFIRPM.

%   Authors: L. Karam, J. McClellan
%   Revised: October 1996, D. Orofino
%
%   Copyright 1988-2004 The MathWorks, Inc.

%  [DH,DW]=HIGHPASS(M,F,GF,W,DELAY)
%       N: filter order (length minus one)
%       F: vector of band edges
%      GF: vector of frequencies at which to evaluate
%       W: vector of weights, one per band
%   DELAY: negative slope of the phase.
%           N/2=(L-1)/2 for exact linear phase.
%
%      DH: vector of desired filter response (mag & phase)
%      DW: vector of weights (positive)
%
% NOTE: DH(GF) and DW(GF) are specified as functions of frequency

% Support query by CFIRPM for the default symmetry option:
if nargin==2,
  % Return symmetry default:
  if strcmp(N,'defaults'),
    % Second arg (F) is cell-array of args passed later to function:
    num_args = length(F);
    % Get the delay value:
    if num_args<5, delay=0; else delay=F{5}; end
    % Use delay arg to base symmetry decision:
    if isequal(delay,0), DH = 'even'; else DH='real'; end
    return
  end
end

% Standard call:
error(nargchk(4,5,nargin,'struct'));
if nargin<5, delay = 0; end
delay = delay + N/2;  % adjust for linear phase

Le = length(F);
if (Le == 4),
  if any(F < 0),
    error(message('signal:highpass:MustBePositive')); 
  end
elseif (Le == 6),
  if F(3)*F(4) > 0,
    error(message('signal:highpass:InvalidFreqVec'));
  end
else
  error(message('signal:highpass:InvalidDimensions'))
end

% Optimization weighting:
W = [1;1]*(W(:).'); W = W(:);

% Construct "highpass" magnitude response:
mags = ones(size(W));
mags(Le-3:Le-2) = 0;  % Zero in 2nd-to-last band

DH = interp1(F(:), mags, GF) .* exp(-1i*pi*GF*delay);
DW = interp1(F(:),    W, GF);

% [EOF]
