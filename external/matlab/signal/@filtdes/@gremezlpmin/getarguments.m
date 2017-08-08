function [F, A, W, args] = getarguments(h, d)
%GETARGUMENTS

%   Author(s): J. Schickler
%   Copyright 1988-2003 The MathWorks, Inc.

[Fpass, Fstop, Dpass, Dstop] = getdesignspecs(h, d);

F = [0 Fpass Fstop 1];
A = [1 1 0 0];
W = [Dpass Dstop];

args = {};

% [EOF]
