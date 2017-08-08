function [order,wn] = buttord(wp,ws,rp,rs,opt)
%BUTTORD Butterworth filter order selection.
%   [N, Wn] = BUTTORD(Wp, Ws, Rp, Rs) returns the order N of the lowest
%   order digital Butterworth filter which has a passband ripple of no more
%   than Rp dB and a stopband attenuation of at least Rs dB. Wp and Ws are
%   the passband and stopband edge frequencies, normalized from 0 to 1
%   (where 1 corresponds to pi radians/sample). For example,
%       Lowpass:    Wp = .1,      Ws = .2
%       Highpass:   Wp = .2,      Ws = .1
%       Bandpass:   Wp = [.2 .7], Ws = [.1 .8]
%       Bandstop:   Wp = [.1 .8], Ws = [.2 .7]
%   BUTTORD also returns Wn, the Butterworth natural frequency (or,
%   the "3 dB frequency") to use with BUTTER to achieve the specifications.
%
%   [N, Wn] = BUTTORD(Wp, Ws, Rp, Rs, 's') does the computation for an
%   analog filter, in which case Wp and Ws are in radians/second.
%
%   When Rp is chosen as 3 dB, the Wn in BUTTER is equal to Wp in BUTTORD.
%
%   % Example 1:
%   %   For data sampled at 1000 Hz, design a lowpass filter with less than
%   %   3 dB of ripple in the passband, defined from 0 to 40 Hz, and at
%   %   least 60 dB of attenuation in the stopband, defined from 150 Hz
%   %   to the Nyquist frequency (500 Hz).
%
%   Wp = 40/500; Ws = 150/500;
%   [n,Wn] = buttord(Wp,Ws,3,60);   % Gives mimimum order of filter
%   [b,a] = butter(n,Wn);           % Butterworth filter design
%   freqz(b,a,512,1000);            % Plots the frequency response
%
%   % Example 2:
%   %   Design a bandpass filter with passband of 60 Hz to 200 Hz, with
%   %   less than 3 dB of ripple in the passband, and 40 dB attenuation in
%   %   the stopbands that are 50 Hz wide on both sides of the passband.
%
%   Wp = [60 200]/500; Ws = [50 250]/500;   % Normalizing frequency
%   Rp = 3; Rs = 40;
%   [n,Wn] = buttord(Wp,Ws,Rp,Rs);  % Gives mimimum order of filter
%   [b,a] = butter(n,Wn);           % Butterworth filter design
%   freqz(b,a,128,1000)             % Plots the frequency response
%
%   See also BUTTER, CHEB1ORD, CHEB2ORD, ELLIPORD, DESIGNFILT.

%   Author(s): L. Shure, 6-9-88
%              T. Krauss, 11-13-92, revised
%   Copyright 1988-2013 The MathWorks, Inc.

%   Reference(s):
%     [1] Rabiner and Gold, p 241.

narginchk(4,5);  
nargoutchk(0,2);

% Cast to enforce precision rules
wp = signal.internal.sigcasttofloat(wp,'double','buttord','Wp',...
  'allownumeric');
ws = signal.internal.sigcasttofloat(ws,'double','buttord','Ws',...
  'allownumeric');
rp = signal.internal.sigcasttofloat(rp,'double','buttord','Rp',...
  'allownumeric');
rs = signal.internal.sigcasttofloat(rs,'double','buttord','Rs',...
  'allownumeric');

if nargin == 4
    opt = 'z';
elseif nargin == 5
    if ~strcmp(opt,'z') && ~strcmp(opt,'s')
        error(message('signal:buttord:InvalidParam'));
    end
end

[msg,msgobj]=freqchk(wp,ws,opt);
if ~isempty(msg), error(msgobj); end

ftype = 2*(length(wp) - 1);
if wp(1) < ws(1)
    ftype = ftype + 1;	% low (1) or reject (3)
else
    ftype = ftype + 2;	% high (2) or pass (4)
end

% first, prewarp frequencies from digital (unit circle) to analog (imag. axis):
if strcmp(opt,'z')	% digital
    WP=tan(pi*wp/2);
    WS=tan(pi*ws/2);
else  % don't have to if analog already
    WP=wp;
    WS=ws;
end
%note - on old systems that are NOT case sensitive, this will still work OK

% next, transform to low pass prototype with passband edge of 1 and stopband
% edges determined by the following: (see Rabiner and Gold, p.258)
if ftype == 1	% low
    WA=WS/WP;
elseif ftype == 2	% high
    WA=WP/WS;
elseif ftype == 3	% stop
    fo = optimset('display','none');
    wp1 = lclfminbnd('bscost',WP(1),WS(1)-1e-12,fo,1,WP,WS,rs,rp,'butter');
    WP(1) = wp1;
    wp2 = lclfminbnd('bscost',WS(2)+1e-12,WP(2),fo,2,WP,WS,rs,rp,'butter');
    WP(2) = wp2;
    WA=(WS*(WP(1)-WP(2)))./(WS.^2 - WP(1)*WP(2));
elseif ftype == 4	% pass
    WA=(WS.^2 - WP(1)*WP(2))./(WS*(WP(1)-WP(2)));
end


% find the minimum order b'worth filter to meet the more demanding spec:
WA=min(abs(WA));
order = ceil( log10( (10 .^ (0.1*abs(rs)) - 1)./ ...
    (10 .^ (0.1*abs(rp)) - 1) ) / (2*log10(WA)) );

% next find the butterworth natural frequency W0 (or, the "3dB frequency")
% to give exactly rs dB at WA.  W0 will be between 1 and WA:
W0 = WA / ( (10^(.1*abs(rs)) - 1)^(1/(2*(abs(order)))));

% now convert this frequency back from lowpass prototype
% to the original analog filter:
if ftype == 1	% low
    WN=W0*WP;
elseif ftype == 2	% high
    WN=WP/W0;
elseif ftype == 3	% stop
    WN(1) = ( (WP(2)-WP(1)) + sqrt((WP(2)-WP(1))^2 + ...
        4*W0.^2*WP(1)*WP(2)))./(2*W0);
    WN(2) = ( (WP(2)-WP(1)) - sqrt((WP(2)-WP(1))^2 + ...
        4*W0.^2*WP(1)*WP(2)))./(2*W0);
    WN=sort(abs(WN));
elseif ftype == 4	% pass
    W0=[-W0 W0];  % need both left and right 3dB frequencies
    WN= -W0*(WP(2)-WP(1))/2 + sqrt( W0.^2/4*(WP(2)-WP(1))^2 + WP(1)*WP(2) );
    WN=sort(abs(WN));
end

% finally, transform frequencies from analog to digital if necessary:
if strcmp(opt,'z')	% digital
    wn=(2/pi)*atan(WN);  % bilinear transform
else
    wn=WN;
end

