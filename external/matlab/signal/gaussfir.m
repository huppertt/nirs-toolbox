function h = gaussfir(BT,NT,OF)
%GAUSSFIR   Gaussian FIR Pulse-Shaping Filter Design.
%
%   WARNING: GAUSSFIR is not recommended. Use GAUSSDESIGN instead.
%
%   H=GAUSSFIR(BT) designs a low pass FIR gaussian pulse-shaping filter.
%   BT is the 3-dB bandwidth-symbol time product where B is the one-sided
%   bandwidth in Hertz and T is in seconds.
%
%   H=GAUSSFIR(BT,NT) NT is the number of symbol periods between the start
%   of the filter impulse response and its peak. If NT is not specified, 
%   NT = 3 is used.
%
%   H=GAUSSFIR(BT,NT,OF) OF is the oversampling factor, that is, the number
%   of samples per symbol. If OF is not specified, OF = 2 is used.
%
%   The length of the impulse response of the filter is given by 2*OF*NT+1.
%   Also, the coefficients H are normalized so that the nominal passband
%   gain is always equal to one.
%
%   % EXAMPLE: Design a Gaussian filter to be used in a GSM GMSK scheme.
%   BT = .3; % 3-dB bandwidth-symbol time
%   OF = 8;  % Oversampling factor (i.e., number of samples per symbol)
%   NT = 2;  % 2 symbol periods to the filters peak. 
%   h = gaussfir(BT,NT,OF); 
%   hfvt = fvtool(h,'impulse');
%
%   See also RCOSDESIGN.

%   References:
%   [1] Rappaport T.S., Wireless Communications Principles and Practice,  
%   2nd Ed., Prentice Hall, 2002.
%   [2] Krishnapura N., Pavan S., Mathiazhagan C., Ramamurthi B., "A
%   Baseband Pulse Shaping Filter for Gaussian Minimum Shift Keying,"
%   Proceedings of the 1998 IEEE International Symposium on Circuits and
%   Systems, 1998. ISCAS '98. 

%   Copyright 2004-2013 The MathWorks, Inc.

% Validate number I/O arguments.
narginchk(1,3);
nargoutchk(0,1);

if nargin < 2, NT = 3; end
if nargin < 3, OF = 2; end

% Check for valid BT
chkBT(BT);
% Cast to enforce precision rules
BT = double(BT);

% Convert to t in which to compute the filter coefficients
t= convert2t(OF,NT);

% Equation 6.53 of [1], page 290 is
% a = sqrt(log(2)/2)/B, here we use alpha = a/T
alpha = sqrt(log(2)/2)/(BT);

% Equation 5.54 of [1] is
% h = (sqrt(pi)/a)*exp(-(t1*pi/a).^2); 
% We use t = t1/T, alpha = a/T.  Then
% h = (sqrt(pi)*T/alpha)*exp(-(t*pi/alpha).^2); 
% But then we normalize, so T is not needed.
h = (sqrt(pi)/alpha)*exp(-(t*pi/alpha).^2); 
 
% Normalize coefficients
h = h./sum(h);


%--------------------------------------------------------------------------
function t = convert2t(OF,NT)

% Check for valid OF and NT
if ~chkInput(OF)
  error(message('signal:gaussfir:invalidOSFactor'));
end
% Cast to enforce precision rules
OF = double(OF);

if ~chkInput(2*NT)
  error(message('signal:gaussfir:invalidNumSPeriods'));
end
% Cast to enforce precision rules
NT = double(NT);

% Filter Length
filtLen = 2*OF*NT+1;
t = linspace(-NT,NT,filtLen);


%-----------------------------------------------------------------------
function chkBT(val)

if isempty(val) || length(val) > 1 || ~isa(val,'numeric') || ...
        ~isreal(val) || val<=0,
    error(message('signal:gaussfir:invalidBT'));
end

%-----------------------------------------------------------------------
function ok = chkInput(val)

ok = true;
if isempty(val) || length(val) > 1 || ~isa(val,'numeric') || ...
        ~isreal(val) || val~=round(val) || val<=0,
      ok = false;
end


% [EOF]
