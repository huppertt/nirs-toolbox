function [DH,DW] = allpass(N, F, GF, W)
%ALLPASS Desired frequency response for allpass filters.
%   CFIRPM(N,F,@allpass, ...) designs a nonlinear phase allpass filter
%   using CFIRPM by adding a quadratic component (4/pi)*(omega)^2 to the
%   usual linear phase found in a symmetric FIR filter .
%
%   See also CFIRPM.

%   Authors: J. McClellan
%   Copyright 1988-2004 The MathWorks, Inc.

error(nargchk(4,4,nargin,'struct'));

DH = exp(-1i*pi*GF*N/2 + 1i*pi*pi*sign(GF).*GF.*GF*(4/pi));
DW = ones(size(GF));

% [EOF] allpass.m
