function f = thisissos(Hd)
%THISISSOS  True if second-order-section.
%   THISISSOS(Hd) returns 1 if filter Hd is second-order or less, and 0
%   otherwise. 
%
%   See also DFILT.   
  
%   Author: Thomas A. Bryan
%   Copyright 1988-2004 The MathWorks, Inc.

% This should be private

f = order(Hd)<=2;
