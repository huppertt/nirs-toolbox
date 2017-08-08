function s = norm(Hd,pnorm,tol)
%NORM   Filter norm.
%   NORM(H) returns the L2-norm of a digital filter (DFILT) H.  
%
%   NORM(H,PNORM) returns the p-norm of a filter. PNORM can be either
%   frequency-domain norms: 'L1', 'L2', 'Linf' or discrete-time-domain
%   norms: 'l1', 'l2', 'linf'. Note that the L2-norm of a filter is equal
%   to its l2-norm (Parseval's theorem), but this is not true for other
%   norms.
%
%   When computing the l1-, l2-, linf-, L1-, and L2-norms of an IIR filter,
%   NORM(...,TOL) will specify the tolerance for greater or less accuracy.
%   By default, TOL = 1e-8.
%
%      EXAMPLE(S):
%      % Compute the L2-norm with a tolerance of 1e-10 for an IIR filter
%      Hs = fdesign.lowpass; % Create a filter design specifications object
%      Hd = butter(Hs);      % Design a Butterworth SOS filter
%      L2 = norm(Hd,'L2',1e-10);% Compute the L2-norm
%
%   See also MFILT/NORM, ADAPTFILT/NORM, DFILT/SCALE, DFILT/SCALECHECK,
%   NORM.

%   Author(s): R. Losada
%   Copyright 2003 The MathWorks, Inc.

error(nargchk(1,3,nargin,'struct'));
if nargin < 2, pnorm = 'L2'; end

% Default tolerance
if nargin < 3, tol = 1e-8; end

% Check for stability
if ~isstable(Hd),
    error(message('signal:dfilt:basefilter:norm:unstableFilter'));
end

switch pnorm,
    case 'Linf',
        % inf-norm is simply given by the max of the magnitude response
        H = freqz(Hd);
        s = max(abs(H));
    case {'L2','l1','l2','linf'}
        s = computetimenorm(Hd,pnorm,tol);
    case 'L1',
        % Inline a function handle used to evaluate the abs of the freq
        % resp.
        f = @(w) abs(freqz(Hd,w));
        s = 1/(2*pi)*quad(f,-pi,pi,tol);
        % Faster, alternate computation using simple rectangle approx to integral
        %Nfft = 8192;
        %[H,w]=freqz(Hd,Nfft,'whole');
        %s = 1/(2*pi)*sum(abs(H))*(2*pi/Nfft);
end

%--------------------------------------------------------------------------
function s = computetimenorm(Hd,pnorm,tol)
%COMPUTETIMENORM   Compute l1-, l2-, and linf-norm of a filter.
%

switch pnorm
    case 'l1',
        p = 1;
    case {'l2','L2'}
        p = 2;
    case {'linf'}
        p = inf;
end

N = impzlength(Hd,tol);
H = impz(Hd,N);
s = norm(H,p);

% [EOF]
