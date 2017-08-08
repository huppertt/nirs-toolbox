function f = thisisscalar(Hd)
%THISISSCALAR  True if scalar filter.
%   THISISSCALAR(Hd) returns 1 if Hd is a scalar filter, and 0 otherwise.
%
%   See also DFILT.   
  
%   Author: Thomas A. Bryan
%   Copyright 1988-2004 The MathWorks, Inc.

% This should be private

f = order(Hd)==0;