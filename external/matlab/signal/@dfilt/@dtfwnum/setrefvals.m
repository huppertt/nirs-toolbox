function setrefvals(this, refvals)
%SETREFVALS   Set reference values.

%   Author(s): R. Losada
%   Copyright 2003 The MathWorks, Inc.

rcnames = refcoefficientnames(this);

set(this,rcnames,refvals);


% [EOF]
