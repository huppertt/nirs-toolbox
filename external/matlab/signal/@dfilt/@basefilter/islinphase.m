function f = islinphase(Hb, varargin)
%ISLINPHASE  True for linear phase filter.
%   ISLINPHASE(Hb) returns 1 if filter Hb is linear phase, and 0 otherwise.
%
%   ISLINPHASE(Hb,TOL) uses tolerance TOL to determine when two numbers are
%   close enough to be considered equal.
%
%   See also DFILT.

%   Authors: J. Schickler
%   Copyright 1988-2002 The MathWorks, Inc.

error(nargchk(1,2,nargin,'struct'));

f = base_is(Hb, 'thisislinphase', varargin{:});

% [EOF]
