function setrefvals(this, refvals)
%SETREFVALS   Set the refvals.

%   Author(s): J. Schickler
%   Copyright 1988-2004 The MathWorks, Inc.

rcnames = refcoefficientnames(this);

set(this,rcnames,refvals);

% [EOF]
