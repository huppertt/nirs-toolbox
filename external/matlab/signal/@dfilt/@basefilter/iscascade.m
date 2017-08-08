function f = iscascade(Hb)
%ISCASCADE  True for cascaded filter.
%   ISCASCADE(Hb) returns 1 if filter Hb is cascade of filters, and 0 otherwise.
%
%   See also DFILT.   

%   Author: J. Schickler
%   Copyright 1988-2002 The MathWorks, Inc.

f = base_is(Hb, 'thisiscascade');

% [EOF]
