function [DH,DW] = multiband(N, F, GF, W, mags, delay)
%MULTIBAND Desired frequency response for multiband filters.
%  CFIRPM(N,F,{@multiband, M} ...) designs a linear-phase multiband filter
%  response using CFIRPM, where M is a vector of magnitudes, one per band
%  edge given in F.
%
%  CFIRPM(N,F,{@multiband, M, D}, ...) specifies group-delay offset D such
%  that the filter response will have a group delay of N/2 + D in units of
%  the sample interval, where N is the filter order. Negative values create
%  less delay, while positive values create more delay.  By default, D=0.
%
%  The symmetry option SYM defaults to 'even' if unspecified in the call to
%  CFIRPM, if no negative band edge frequencies are specified in F.
%
%  EXAMPLE: Design a 31-tap, complex multiband filter.
%    b = cfirpm(30,[-1 -.5 -.4 .7 .8 1],{@multiband,[0 0 1 2 0 0]});
%    fvtool(b);  % View filter response.
%
%  EXAMPLE: Reduced group delay filter response.
%    b = cfirpm(30,[-1 -.5 -.4 .7 .8 1],{@multiband,[0 0 1 2 0 0],-1});
%    fvtool(b);
%
%  See also CFIRPM.

%   Authors: L. Karam, J. McClellan
%   Revised: October 1996, D. Orofino
%
%   Copyright 1988-2004 The MathWorks, Inc.

%  [DH,DW]=MULTIBAND(N,F,GF,W,MAGS,DELAY)
%       N: filter order (length minus one)
%       F: vector of band edges
%      GF: vector of frequencies at which to evaluate
%       W: vector of weights, one per band
%    MAGS: vector of mags, one per band edge
%   DELAY: negative slope of the phase.
%           N/2=(L-1)/2 for exact linear phase.
%
%     DH: vector of desired filter response (mag & phase)
%     DW: vector of weights (positive)
%
% NOTE: DH(GF) and DW(GF) are specified as functions of frequency

% Ex. 1: Asymmetric complex response:
%   b=cfirpm(32,[-.8 -.5 -.3 .2 .4 .6],[1 1 0 0 1 1]);
%

% Support query by CFIRPM for the default symmetry option:
if nargin==2,
  % Return symmetry default:
  if strcmp(N,'defaults'),
    % Second arg (F) is cell-array of args passed later to function:
    num_args = length(F);
    % Get the delay value:
    if num_args<6, delay=0; else delay=F{6}; end
    % Use delay arg to base symmetry decision:
    if isequal(delay,0), DH='even'; else DH='real'; end
    return
  end
end

% Standard call:
error(nargchk(4,6,nargin,'struct'));
if nargin<5,
  error(message('signal:multiband:InvalidParam'));
end
if nargin<6, delay = 0; end
delay = delay + N/2;  % adjust for linear phase

W  = [1;1] * (W(:).');
DH = interp1(F(:), mags(:), GF) .* exp(-1i*pi*GF*delay);
DW = interp1(F(:),    W(:), GF);

% [EOF]
