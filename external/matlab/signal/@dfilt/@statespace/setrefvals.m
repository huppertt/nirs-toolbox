function setrefvals(this, refvals)
%SETREFVALS   Set reference values.
%This should be a private method.

%   Author(s): R. Losada
%   Copyright 2003-2004 The MathWorks, Inc.

rcnames = refcoefficientnames(this);

set(this,rcnames,refvals);


% [EOF]
