function f = issos(Hb)
%ISSOS  True if second-order-section.
%   ISSOS(Hb) returns 1 if filter Hb is second-order or less, and 0
%   otherwise. 
%
%   See also DFILT.   
  
%   Author: J. Schickler
%   Copyright 1988-2002 The MathWorks, Inc.

f = base_is(Hb, 'thisissos');

% [EOF]
