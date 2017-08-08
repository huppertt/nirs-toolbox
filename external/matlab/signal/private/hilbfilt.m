function [DH,DW] = hilbfilt(N, F, GF, W, delay)
%HILBFILT Desired frequency response for Hilbert filters.
%  CFIRPM(N,F,@hilbfilt, ...) designs a linear-phase hilbert filter
%  response using CFIRPM.
%
%  CFIRPM(N,F,{@hilbfilt, D}, ...) specifies group-delay offset D such
%  that the filter response will have a group delay of N/2 + D in
%  units of the sample interval, where N is the filter order.
%  Negative values create less delay, while positive values create
%  more delay.  By default, D=0.
%
%  Note that DC must be in a transition band.
%
%  The symmetry option SYM defaults to 'odd' if unspecified in the
%  call to CFIRPM, if no negative band edge frequencies are
%  specified in F.
%
%  EXAMPLE: Determine analytic envelope of a signal
%     N   = 32; t = (0:200)';
%     env = sin(pi*t/200);          % Signal envelope
%     x   = sin(2*pi*0.3*t) .* env; % Test signal
%     b   = cfirpm(32,[.1 .9],@hilbfilt);
%     pad = zeros(N/2,1);           % Delay compensation
%     y   = filter(b,1,[x;pad]);    % Hilbert transform
%     as  = [pad;x] + 1i*y;         % Analytic signal
%     ae  = abs(as(N/2+1:end));     % Analytic envelope
%     plot(t,ae, t,env);
%
%  See also CFIRPM.

%   Authors: L. Karam, J. McClellan
%   Revised: October 1996, D. Orofino
%
%   Copyright 1988-2010 The MathWorks, Inc.

%  [DH,DW]=HILBFILT(N,F,GF,W,DELAY)
%       N: filter order (length minus one)
%       F: vector of band edges
%      GF: vector of frequencies at which to evaluate
%       W: vector of weights, one per band
%   DELAY: negative slope of the phase.
%           N/2=(L-1)/2 for exact linear phase.
%
%     DH: vector of desired filter response (mag & phase)
%     DW: vector of weights (positive)
%
% NOTE: DH(f) and DW(f) are specified as functions of frequency

% Support query by CFIRPM for the default symmetry option:
if nargin==2,
  % Return symmetry default:
  if strcmp(N,'defaults'),
    % Second arg (F) is cell-array of args passed later to function:
    num_args = length(F);
    % Get the delay value:
    if num_args<5, delay=0; else delay=F{5}; end
    % Use delay arg to base symmetry decision:
    if isequal(delay,0), DH='odd'; else DH='real'; end
    return
  end
end

% Standard call:
error(nargchk(4,5,nargin,'struct'));
if nargin < 5, delay = 0; end
delay = delay + N/2;  % adjust for linear phase

Le = length(F);
if Le == 2,
  if any(F <= 0),
    error(message('signal:hilbfilt:MustBePositive'));
  end
elseif Le == 4,
  if F(2)*F(3) > 0,
    error(message('signal:hilbfilt:InvalidParam'));
  end
else
  error(message('signal:hilbfilt:InvalidDimensions'))
end

W    = [1;1] * (W(:).'); W = W(:);
mags = ones(size(W));
DH   = interp1(F(:), mags, GF) .* -1i .* exp(-1i*pi*GF*delay);
DW   = interp1(F(:),    W, GF);

% end of hilbfilt.m
