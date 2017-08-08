function f = thisisreal(Hd)
%ISREAL  True for filter with real coefficients.
%   ISREAL(Hd) returns 1 if filter Hd has real coefficients, and 0
%   otherwise. 
%
%   See also DFILT.   
  
%   Author: Thomas A. Bryan, J. Schickler
%   Copyright 1988-2004 The MathWorks, Inc.

% This should be private

f = all(isreal(Hd.Stage));

