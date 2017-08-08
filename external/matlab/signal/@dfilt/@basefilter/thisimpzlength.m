function n = thisimpzlength(this, varargin)
%THISIMPZLENGTH   Dispatch and call the method.

%   Author(s): J. Schickler
%   Copyright 1988-2004 The MathWorks, Inc.

Hd = dispatch(this);
n = impzlength(Hd, varargin{:});

% [EOF]
