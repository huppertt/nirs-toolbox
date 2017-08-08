function h = alpmin(wp,ws,rp,rs)
%ALPMIN   Construct an ALPMIN object.
%   H = ALPMIN(Wpass,Wstop,Apass,Astop) constructs an analog minimum-order
%   lowpass filter specifications object.
%
%   Wpass is the passband-edge frequency in radians-per-second and must be
%   a positive scalar.
%
%   Wstop is the stopband-edge frequency in radians-per-second and must be
%   a positive scalar.
%
%   Apass is the maximum passband deviation in dB. It must be a positive
%   scalar.
%
%   Astop is the minimum stopband attenuation in dB. It must be a positive
%   scalar.

%   Author(s): R. Losada
%   Copyright 1988-2003 The MathWorks, Inc.

h = fspecs.alpmin;

h.ResponseType = 'Analog minimum-order lowpass';
if nargin > 0,
    h.Wpass = wp;
end

if nargin > 1,
    h.Wstop = ws;
end

if nargin > 2,
    h.Apass = rp;
end

if nargin > 3,
    h.Astop = rs;
end


% [EOF]
