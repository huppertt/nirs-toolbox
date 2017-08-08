function s = saveobj(this)
%SAVEOBJ  Save this object.
%   OUT = SAVEOBJ(ARGS) <long description>

%   Copyright 2006 The MathWorks, Inc.

s.class   = class(this);

% Save all of the public properties.
s = setstructfields(s, get(this));

% [EOF]
