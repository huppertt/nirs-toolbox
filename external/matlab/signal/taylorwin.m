function w = taylorwin(varargin)
%TAYLORWIN   Taylor window.
%
% TAYLORWIN(N) returns an N-point Taylor window in a column vector.
%
% TAYLORWIN(N,NBAR) returns an N-point Taylor window with NBAR nearly
% constant-level sidelobes adjacent to the mainlobe. NBAR must be an
% integer greater than or equal to one.
%
% TAYLORWIN(N,NBAR,SLL) returns an N-point Taylor window with SLL maximum
% sidelobe level in dB relative to the mainlobe peak. SLL must be a
% negative value, e.g., -30 dB.
%
% NBAR should satisfy NBAR >= 2*A^2+0.5, where A is equal to
% acosh(10^(-SLL/20))/pi, otherwise the sidelobe level specified is not
% guaranteed. If NBAR is not specified it defaults to 4. SLL defaults to
% -30 dB if not specified.
%
% EXAMPLE
%   This example generates a 64-point Taylor window with 4 sidelobes
%   adjacent to the mainlobe that are nearly constant-level, and a peak
%   sidelobe level of -35 dB relative to the mainlobe peak.
%
%   w = taylorwin(64,5,-35);
%   wvtool(w);
%
%   See also CHEBWIN.

%   Copyright 2005-2013 The MathWorks, Inc.

%
%   References:
%     [1] Carrara, Walter G., Ronald M. Majewski, and Ron S. Goodman,
%         Spotlight Synthetic Aperture Radar: Signal Processing Algorithms,
%         Artech House, October 1, 1995.
%     [2] Brookner, Eli, Practical Phased Array Antenna Systems,
%         Artech House, Inc., 1991, pg. 2-51.

% Validate input and set default values.

narginchk(1, 3);  
[N,NBAR,SLL] = validateinputs(varargin{:});

A = acosh(10^(-SLL/20))/pi;

% Taylor pulse widening (dilation) factor.
sp2 = NBAR^2/(A^2 + (NBAR-.5)^2);

w = ones(N,1);
Fm = zeros(NBAR-1,1);
summation = 0;
k = (0:N-1)';
xi = (k-0.5*N+0.5)/N;
for m = 1:(NBAR-1),
    Fm(m) = calculateFm(m,sp2,A,NBAR);
    summation = Fm(m)*cos(2*pi*m*xi)+summation;
end
w = w + 2*summation;

%-------------------------------------------------------------------
function Fm = calculateFm(m,sp2,A,NBAR)
% Calculate the cosine weights.

n = (1:NBAR-1)';
p = [1:m-1, m+1:NBAR-1]'; % p~=m

Num = prod((1 - (m^2/sp2)./(A^2+(n-0.5).^2)));
Den = prod((1 - m^2./p.^2));

Fm = ((-1)^(m+1).*Num)./(2.*Den);

%-------------------------------------------------------------------
function [N,NBAR,SLL] = validateinputs(varargin)
% Cast to enforce Precision Rules
N = signal.internal.sigcasttofloat(varargin{1},'double','taylorwin',...
  'N','allownumeric');
% Validate order
if N ~= floor(N),
    N = round(N);
    warning(message('signal:taylorwin:WindowLengthMustBeInteger'));
end

if  N < 0,
    error(message('signal:taylorwin:WindowLengthMustBePositive'));
end

% Validate NBAR
if nargin < 2,
    NBAR = 4;
else
    % Cast to enforce Precision Rules
    NBAR = signal.internal.sigcasttofloat(varargin{2},'double',...
      'taylorwin','NBAR','allownumeric');
    if NBAR<=0 || NBAR~=floor(NBAR),
        error(message('signal:taylorwin:NBARMustBeInteger'));
    end
end

% Validate SLL
if nargin < 3,
    SLL = -30;
else
    % Cast to enforce Precision Rules
    SLL = signal.internal.sigcasttofloat(varargin{3},'double','taylorwin',...
      'SLL','allownumeric');    
    if SLL>0,
        error(message('signal:taylorwin:SLLMustBeNegative'));
    end
end

% [EOF]
