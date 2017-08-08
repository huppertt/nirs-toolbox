function f = isreal(Hb)
%ISREAL  True for filter with real coefficients.
%   ISREAL(Hb) returns 1 if filter Hb has real coefficients, and 0
%   otherwise. 
%
%   See also DFILT.   
  
%   Author: J. Schickler
%   Copyright 1988-2002 The MathWorks, Inc.

% Use strings because construct of fcn handles to functions not on the path
% (methods) is very slow
f = base_is(Hb, 'thisisreal');

% [EOF]
