function f = thisissos(Hd)
%ISSOS  True if second-order-section.
%   ISSOS(Hd) returns 1 if filter Hd is second-order or less, and 0
%   otherwise. 
%
%   See also DFILT.   
  
%   Author: Thomas A. Bryan, J. Schickler
%   Copyright 1988-2004 The MathWorks, Inc.

% This should be private

f = all(issos(Hd.Stage));

