function h = alppass(N,Wp,Apass)
%ALPPASS   Construct an ALPPASS object.
%   H = ALPPASS(N,Wp,Apass) constructs an analog lowpass filter design
%   specifications object H with passband-edge specs.
%
%   N is the filter order, and must be a positive integer.
%
%   Wp is the passband-edge frequency, in radians-per-second.
%
%   Apass is the maximum passband ripple, in dB.

%   Author(s): R. Losada
%   Copyright 1988-2003 The MathWorks, Inc.

h = fspecs.alppass;

h.ResponseType = 'Analog lowpass with passband-edge specifications';
if nargin > 0,
    h.FilterOrder = N;
end

if nargin > 1,
    h.Wpass = Wp;
end

if nargin > 2,
    h.Apass = Apass;
end

% [EOF]
