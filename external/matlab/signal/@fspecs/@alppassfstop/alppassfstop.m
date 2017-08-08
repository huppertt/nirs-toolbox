function h = alppassfstop(N,Wp,Ws,Apass)
%ALPPASSFSTOP   Construct an ALPPASSFSTOP object.
%   H = ALPPASSFSTOP(N,Wp,Ws,Apass) constructs an analog lowpass filter
%   design specifications object H with passband-edge specs and
%   stopband-edge frequency.
%
%   N is the filter order, and must be a positive integer.
%
%   Wp is the passband-edge frequency, in radians-per-second.
%
%   Ws is the stopband-edge frequency, in radians-per-second. It must be
%   larger than Wp.
%
%   Apass is the maximum passband ripple, in dB.

%   Author(s): R. Losada
%   Copyright 1988-2003 The MathWorks, Inc.

h = fspecs.alppassfstop;

h.ResponseType = 'Analog lowpass with passband-edge specifications and stopband frequency';
if nargin > 0,
    h.FilterOrder = N;
end

if nargin > 1,
    h.Wpass = Wp;
end

if nargin > 2,
    h.Wstop = Ws;
end

if nargin > 3,
    h.Apass = Apass;
end

% [EOF]
