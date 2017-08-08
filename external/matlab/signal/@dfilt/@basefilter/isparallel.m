function f = isparallel(Hb)
%ISPARALLEL  True for filter with parallel sections.
%   ISPARALLEL(Hb) returns 1 if filter Hb is composed of parallel sections,
%   and 0 otherwise. 
%
%   See also DFILT.   
  
%   Author: J. Schickler
%   Copyright 1988-2002 The MathWorks, Inc.

f = base_is(Hb, 'thisisparallel');

% [EOF]
