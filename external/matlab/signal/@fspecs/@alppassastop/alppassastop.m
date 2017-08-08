function h = alppassastop(N,Wp,Apass,Astop)
%ALPPASSASTOP   Construct an ALPPASSASTOP object.
%   H = ALPPASSASTOP(N,Wp,Apass,Astop) constructs an analog lowpass filter
%   design specifications object H with passband-edge specs and stopband-edge
%   frequency.
%
%   N is the filter order, and must be a positive integer.
%
%   Wp is the passband-edge frequency, in radians-per-second.
%
%   Apass is the maximum passband ripple, in dB.
%
%   Astop is the minimum stopband attenuation, in dB.


%   Author(s): R. Losada
%   Copyright 1988-2003 The MathWorks, Inc.

h = fspecs.alppassastop;

h.ResponseType = 'Analog lowpass with passband-edge specifications and stopband attenuation';
if nargin > 0,
    h.FilterOrder = N;
end

if nargin > 1,
    h.Wpass = Wp;
end

if nargin > 2,
    h.Apass = Apass;
end

if nargin > 3,
    h.Astop = Astop;
end

% [EOF]
