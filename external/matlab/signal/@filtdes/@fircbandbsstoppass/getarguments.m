function [F, A, W, args] = getarguments(h, d)
%GETARGUMENTS Return the design method arguments

%   Author(s): J. Schickler
%   Copyright 1988-2003 The MathWorks, Inc.

F = [0 get(d, 'Fpass1') get(d, 'Fstop1') get(d, 'Fstop2') get(d, 'Fpass2') 1];
A = [1 1 0 0 1 1];

mu = get(d, 'magUnits'); set(d, 'magUnits', 'linear');
W = [1 get(d, 'Dstop') get(d, 'Dpass2')]; set(d, 'magUnits', mu);

args = {};

% [EOF]
