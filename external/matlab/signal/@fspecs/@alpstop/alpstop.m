function h = alpstop(N,Ws,Astop)
%ALPSTOP   Construct an ALPSTOP object.
%   H = ALPSTOP(N,Ws,Astop) constructs an analog lowpass filter design
%   specifications object H with stopband-edge specs.
%
%   N is the filter order, and must be a positive integer.
%
%   Ws is the stopband-edge frequency, in radians-per-second.
%
%   Astop is the minimum stopband attenuation, in dB.

%   Author(s): R. Losada
%   Copyright 1988-2003 The MathWorks, Inc.

h = fspecs.alpstop;

h.ResponseType = 'Analog lowpass with stopband-edge specifications';
if nargin > 0,
    h.FilterOrder = N;
end

if nargin > 1,
    h.Wstop = Ws;
end

if nargin > 2,
    h.Astop = Astop;
end

% [EOF]
