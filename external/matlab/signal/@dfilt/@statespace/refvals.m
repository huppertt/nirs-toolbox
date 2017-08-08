function rcvals = refvals(this)
%REFVALS   Return the reference values.

%   Author(s): J. Schickler
%   Copyright 1988-2004 The MathWorks, Inc.

rcnames = refcoefficientnames(this);

rcvals = get(this,rcnames);

% [EOF]
