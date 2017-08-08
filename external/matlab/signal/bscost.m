function n = bscost(wp,ind,WP,WS,rs,rp,filtype)
%BSCOST  Band Stop Cost function for order minimization w.r.t passband edge.
%   BSCOST(wp,ind,WP,WS,rs,rp,FILTYPE) returns the order (non-integer in general)
%   for an analog band stop filter.  This is used by the order estimation
%   routines in minimizing the filter order.
%   Inputs:
%     wp - varying passband edge.
%     ind - index specifying which passband edge is being varied (1 == lower,
%           2 == upper).
%     WP - two element vector of fixed passband edges.
%     WS - two element vector of fixed stopband edges.
%     rs - amount in dB of attenuation in the stopband.
%     rp - amount in dB of ripple in the passband.
%     FILTYPE - can be 'butter', 'cheby', or 'ellip'.
%   Output:
%     n - filter order (fractional).
%
%   See also BUTTORD, CHEB1ORD, CHEB2ORD, ELLIPORD.

%   Author(s): T. Krauss, 2-20-96
%   Copyright 1988-2004 The MathWorks, Inc.

WP(ind) = wp;

% This is old code, so we don't want to remove it, but its immediately
% overwritten by the next line.
% WA=([ WS(1) -WS(2) ] *(WP(1)-WP(2)))./([WS(1) -WS(2)].^2 - WP(1)*WP(2));
WA=(WS*(WP(1)-WP(2)))./(WS.^2 - WP(1)*WP(2));

% find the minimum order filter to meet the more demanding spec:
WA=min(abs(WA));
switch filtype

case 'butter'
    n = ( log10( (10 .^ (0.1*abs(rs)) - 1)./ ...
        (10 .^ (0.1*abs(rp)) - 1) ) / (2*log10(WA)) );
case 'cheby'
    n=(acosh(sqrt((10^(.1*abs(rs))-1)/(10^(.1*abs(rp))-1)))/acosh(WA));
case 'ellip'
    epsilon = sqrt(10^(0.1*rp)-1);
    k1 = epsilon/sqrt(10^(0.1*rs)-1);
    k = 1/WA;
    capk = ellipke([k^2 1-k^2]);
    capk1 = ellipke([(k1^2) 1-(k1^2)]);
    n = (capk(1)*capk1(2)/(capk(2)*capk1(1)));
end

