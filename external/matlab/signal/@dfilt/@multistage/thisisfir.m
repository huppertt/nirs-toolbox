function f = thisisfir(Hd)
%THISISFIR  True for FIR filter.
%   THISISFIR(Hd) returns 1 if filter Hd is FIR, and 0 otherwise.
%
%   See also DFILT.

%   Author: Thomas A. Bryan, J. Schickler
%   Copyright 1988-2004 The MathWorks, Inc.

% This should be private

f = all(isfir(Hd.Stage));
