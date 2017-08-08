function n = thisorder(Hd)
%THISORDER Filter order.
%   THISORDER(Hd) returns the order of filter Hd.
%
%   See also DFILT.   
  
%   Copyright 1988-2012 The MathWorks, Inc.

% This should be private

[b,a] = tf(Hd);
n = filtord(b,a);
