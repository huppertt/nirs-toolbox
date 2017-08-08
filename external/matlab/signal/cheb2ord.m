function [order,wn] = cheb2ord(wp,ws,rp,rs,opt)
%CHEB2ORD Chebyshev Type II filter order selection.
%   [N, Ws] = CHEB2ORD(Wp, Ws, Rp, Rs) returns the order N of the lowest
%   order digital Chebyshev Type II filter which has a passband ripple of
%   no more than Rp dB and a stopband attenuation of at least Rs dB. Wp and
%   Ws are the passband and stopband edge frequencies, normalized from 0 to
%   1 (where 1 corresponds to pi radians/sample). For example,
%       Lowpass:    Wp = .1,      Ws = .2
%       Highpass:   Wp = .2,      Ws = .1
%       Bandpass:   Wp = [.2 .7], Ws = [.1 .8]
%       Bandstop:   Wp = [.1 .8], Ws = [.2 .7]
%   CHEB2ORD also returns Wst, the Chebyshev natural frequency to use with
%   CHEBY2 to achieve the specifications.
%
%   [N, Ws] = CHEB2ORD(Wp, Ws, Rp, Rs, 's') does the computation for an
%   analog filter, in which case Wp and Ws are in radians/second.
%
%   % Example 1:
%   %   For data sampled at 1000 Hz, design a lowpass filter with less than
%   %   3 dB of ripple in the passband defined from 0 to 40 Hz, and at
%   %   least 60 dB of attenuation in the stopband defined from 150 Hz to
%   %   the Nyquist frequency (500 Hz).
%
%   Wp = 40/500; Ws = 150/500;
%   Rp = 3; Rs = 60;
%   [n,Ws] = cheb2ord(Wp,Ws,Rp,Rs)  % Gives mimimum order of filter
%   [b,a] = cheby2(n,Rs,Ws);        % Chebyshev Type II filter
%   freqz(b,a,512,1000);            % Plots the frequency response
%
%   % Example 2:
%   %   Design a bandpass filter with a passband of 60 Hz to 200 Hz, with
%   %   less than 3 dB of ripple in the passband, and 40 dB attenuation in
%   %   the stopbands that are 50 Hz wide on both sides of the passband.
%
%   Wp = [60 200]/500; Ws = [50 250]/500;
%   Rp = 3; Rs = 40;
%   [n,Ws] = cheb2ord(Wp,Ws,Rp,Rs)      % Gives mimimum order of filter
%   [b,a] = cheby2(n,Rs,Ws);            % Chebyshev Type II filter
%   freqz(b,a,512,1000)                 % Plots the frequency response
%
%   See also CHEBY2, CHEB1ORD, BUTTORD, ELLIPORD, DESIGNFILT.

%   Reference: Rabiner and Gold, p 241.

%   Copyright 1988-2013 The MathWorks, Inc.

narginchk(4,5);   
nargoutchk(0,2); 

% Cast to enforce precision rules
wp = signal.internal.sigcasttofloat(wp,'double','cheb2ord','Wp',...
  'allownumeric');
ws = signal.internal.sigcasttofloat(ws,'double','cheb2ord','Ws',...
  'allownumeric');
rp = signal.internal.sigcasttofloat(rp,'double','cheb2ord','Rp',...
  'allownumeric');
rs = signal.internal.sigcasttofloat(rs,'double','cheb2ord','Rs',...
  'allownumeric');

if nargin == 4
    opt = 'z';
elseif nargin == 5
    if ~strcmp(opt,'z') && ~strcmp(opt,'s')
        error(message('signal:cheb2ord:InvalidParam'));
    end
end

[msg,msgobj]=freqchk(wp,ws,opt);
if ~isempty(msg), error(msgobj); end

% figure out filter type
ftype = 2*(length(wp) - 1);
if wp(1) < ws(1)
    ftype = ftype + 1;	% low (1) or reject (3)
else
    ftype = ftype + 2;	% high (2) or pass (4)
end

% first, prewarp frequencies from digital (unit circle) to analog (imag. axis):
if strcmp(opt,'z')	% digital
    WPA=tan(pi*wp/2);
    WSA=tan(pi*ws/2);
else  % don't have to if analog already
    WPA=wp;
    WSA=ws;
end

% next, transform to low pass prototype with passband edge of 1 and stopband
% edges determined by the following: (see Rabiner and Gold, p.258)
if ftype == 1	% low
    WA=WSA/WPA;
elseif ftype == 2	% high
    WA=WPA/WSA;
elseif ftype == 3	% stop
    fo = optimset('display','none');
    wp1 = lclfminbnd('bscost',WPA(1),WSA(1)-1e-12,fo,1,...
        WPA,WSA,rs,rp,'cheby');
    WPA(1) = wp1;
    wp2 = lclfminbnd('bscost',WSA(2)+1e-12,WPA(2),fo,2,...
        WPA,WSA,rs,rp,'cheby');
    WPA(2) = wp2;
    WA=(WSA*(WPA(1)-WPA(2)))./(WSA.^2 - WPA(1)*WPA(2));
elseif ftype == 4	% pass
    WA=(WSA.^2 - WPA(1)*WPA(2))./(WSA*(WPA(1)-WPA(2)));
end

% find the minimum order cheby. type 2 filter to meet the more demanding spec:
WA=min(abs(WA));
order=ceil(acosh(sqrt((10^(.1*abs(rs))-1)/(10^(.1*abs(rp))-1)))/acosh(WA));
% ref: M.E. Van Valkenburg, "Analog Filter Design", p.232, eqn 8.39

% next find the frequency "new_wp" at which the response of an analog low-pass
% prototype is exactly -rp dB.  The prototype is the one for which the beginning
% of the stop-band is at frequency 1.
% (new_wp will be less than 1):
new_wp=1/cosh(1/order*acosh(sqrt((10^(.1*abs(rs)) - 1)/(10^(.1*abs(rp)) - 1))));

% Now convert the stop band frequency back from lowpass prototype
% to the original analog filter.
% Here we use the mapping which maps the frequency "new_wp" to the original WP,
% to map the (+/-)1 frequency to WN (WN will be between WP and WS):
if ftype == 1	% low
    WN=WPA/new_wp;
elseif ftype == 2	% high
    WN=WPA*new_wp;
elseif ftype == 3	% stop
    WN=(WPA(1)-WPA(2))*new_wp/2 + ...
        sqrt( (WPA(2)-WPA(1))^2*new_wp^2/4 + WPA(1)*WPA(2));
    WN(2)=WPA(1)*WPA(2)/WN(1);
elseif ftype == 4	% pass
    WN=(WPA(1)-WPA(2))/(2*new_wp) + ...
        sqrt( (WPA(2)-WPA(1))^2/(4*new_wp^2) + WPA(1)*WPA(2));
    WN(2)=WPA(1)*WPA(2)/WN(1);
    %  WA=(WP.^2 - WN(1)*WN(2))./(WP*(WN(2)-WN(1))) <--- to check, this should be
    %   	-/+ new_wp
end

% finally, transform frequencies from analog to digital if necessary:
if strcmp(opt,'z')	% digital
    wn = ws;  % bilinear transform
else
    wn = WN;
end
