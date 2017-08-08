function f = isscalarstructure(Hb)
%ISSCALARSTRUCTURE  True if scalar filter.
%   ISSCALARSTRUCTURE(Hb) returns 1 if Hb is a scalar filter structure, and 0
%   otherwise. 
%
%   See also DFILT.   
  
%   Author: J. Schickler
%   Copyright 1988-2002 The MathWorks, Inc.

f = base_is(Hb, 'thisisscalarstructure');

% [EOF]
