function [F, A, W, args] = getarguments(h, d)
%GETARGUMENTS Returns the design function arguments

%   Author(s): J. Schickler
%   Copyright 1988-2003 The MathWorks, Inc.

[Fstop1, Fpass1, Fpass2, Fstop2, Dstop1, Dpass, Dstop2] = getdesignspecs(h, d);

F = [0 Fstop1 Fpass1 Fpass2 Fstop2 1];
A = [0 0 1 1 0 0];
W = [Dstop1 Dpass Dstop2];

args = {};

% [EOF]
