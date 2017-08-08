function h = ahpmin(ws,wp,rs,rp)
%AHPMIN   Construct an AHPMIN object.
%   H = AHPMIN(Wstop,Wpass,Astop,Apass) constructs an analog minimum-order
%   lowpass filter specifications object.
%
%   Wstop is the stopband-edge frequency in radians-per-second and must be
%   a positive scalar.
%
%   Wpass is the passband-edge frequency in radians-per-second and must be
%   a positive scalar.
%
%   Astop is the minimum stopband attenuation in dB. It must be a positive
%   scalar.
%
%   Apass is the maximum passband deviation in dB. It must be a positive
%   scalar.



%   Author(s): R. Losada
%   Copyright 1988-2003 The MathWorks, Inc.

h = fspecs.ahpmin;

h.ResponseType = 'Analog minimum-order highpass';

if nargin > 0,
    h.Wstop = ws;
else
    % Override lowpass default
    h.Wstop = 7;
end

if nargin > 1,
    h.Wpass = wp;
else
    % Override lowpass default
    h.Wpass = 13;
end

if nargin > 2,
    h.Astop = rs;
end

if nargin > 3,
    h.Apass = rp;
end



% [EOF]
