function h = gaussdesign(bt,span,sps)
%GAUSSDESIGN Gaussian FIR Pulse-Shaping Filter Design
%   H = GAUSSDESIGN(BT,SPAN,SPS) designs a lowpass FIR Gaussian
%   pulse-shaping filter. BT is the 3-dB bandwidth-symbol time product,
%   where B is the one-sided bandwidth in hertz and T is the symbol time in
%   seconds. The filter is truncated to SPAN symbols and each symbol is
%   represented by SPS samples.
% 
%   GAUSSDESIGN designs a symmetric filter. Therefore, the filter order,
%   which is SPS*SPAN, must be even. Also, the coefficients, H, are
%   normalized so that the nominal passband gain is always equal to one.
%
%   H = GAUSSDESIGN(BT) is the same as GAUSSDESIGN(BT,SPAN,SPS) with SPAN
%   set to 3 and SPS set to 2.
%
%   H = GAUSSDESIGN(BT,SPAN) is the same as GAUSSDESIGN(BT,SPAN,SPS) with
%   SPS set to 2.
% 
%   % EXAMPLE: Design a Gaussian filter to be used in a GSM GMSK scheme.
%   bt   = 0.3; % 3-dB bandwidth-symbol time
%   sps  = 8;   % samples per symbol
%   span = 4;   % filter spans 4 symbols
%   h = gaussdesign(bt,span,sps); 
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

%   Copyright 2013 The MathWorks, Inc.

% Validate number I/O arguments.
narginchk(1,3);
nargoutchk(0,1);

% Check for validity of inputs
validateattributes(bt, {'double','single'}, ...
  {'scalar','nonempty','nonnan','real','positive'}, ...
  'gaussdesign', 'BT', 1);

if nargin >= 2
    validateattributes(span, {'double','single'}, ...
      {'scalar','nonempty','nonnan','real','positive','integer'}, ...
      'gaussdesign', 'SPAN', 2);
else
    span = 3;
end

if nargin == 3
    validateattributes(sps, {'double','single'}, ...
      {'scalar','nonempty','nonnan','real','positive','integer'}, ...
      'gaussdesign', 'SPS', 3);
else
    sps = 2;
end

if mod(sps*span, 2) == 1
  error(message('signal:rcosdesign:OddFilterOrder', sps, span))
end

% Convert to t in which to compute the filter coefficients
filtLen = sps*span+1;
t = linspace(-span/2,span/2,filtLen);

% Equation 6.53 of [1], page 290 is
% a = sqrt(log(2)/2)/B, here we use alpha = a/T
alpha = sqrt(log(2)/2)/(bt);

% Equation 5.54 of [1] is
% h = (sqrt(pi)/a)*exp(-(t1*pi/a).^2); 
% We use t = t1/T, alpha = a/T.  Then
% h = (sqrt(pi)*T/alpha)*exp(-(t*pi/alpha).^2); 
% But then we normalize, so T is not needed.
h = (sqrt(pi)/alpha)*exp(-(t*pi/alpha).^2); 
 
% Normalize coefficients
h = h./sum(h);
