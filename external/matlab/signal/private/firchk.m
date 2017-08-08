function [n,msg1,msg2,msgobj] = firchk(n,Fend,a,exception)
%FIRCHK   Check if specified filter order is valid.
%   FIRCHK(N,Fend,A) checks if the specified order N is valid given the
%   final frequency point Fend and the desired magnitude response vector A.
%   Type 2 linear phase FIR filters (symmetric, odd order) must have a
%   desired magnitude response vector that ends in zero if Fend = 1.  This
%   is because type 2 filters necessarily have a zero at w = pi.
%
%   If the order is not valid, a warning is given and the order
%   of the filter is incremented by one.
%
%   If A is a scalar (as when called from fircls1), A = 0 is
%   interpreted as lowpass and A = 1 is interpreted as highpass.
%
%   FIRCHK(N,Fend,A,EXCEPTION) will not warn or increase the order
%   if EXCEPTION = 1.  Examples of EXCEPTIONS are type 4 filters
%   (such as differentiators or hilbert transformers) or non-linear
%   phase filters (such as minimum and maximum phase filters).

%   Author : R. Losada
%   Copyright 1988-2004 The MathWorks, Inc.

error(nargchk(3,4,nargin,'struct'));

if nargin == 3,
    exception = false;
end

msg1 = '';
msg2 = '';
msgobj = [];
oddord = false; % Flag, initially we assume even order

if isempty(n) || length(n) > 1 || ~isnumeric(n) || ~isreal(n) || n~=round(n) || n<=0,
    msgobj = message('signal:firchk:NeedRealPositiveOrder');
    msg1 = getString(msgobj);
    return
end

if rem(n,2) == 1,
    oddord = true; % Overwrite flag
end
 
if (a(end) ~= 0) && Fend == 1 && oddord && ~exception,
    msgobj = message('signal:firchk:NeedZeroGain');
    msg2 = getString(msgobj);
    n = n+1;
end
    

