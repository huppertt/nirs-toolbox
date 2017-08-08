function b = thisisscalar(this)
%THISISSCALAR   Returns true if all the sections are scalar.

%   Author(s): J. Schickler
%   Copyright 1988-2004 The MathWorks, Inc.

b = all(isscalar(this.Stage));

% [EOF]
