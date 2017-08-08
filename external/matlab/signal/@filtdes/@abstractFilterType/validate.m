function varargout = validate(this, d)
%VALIDATE Returns true if this object is valid

%   Copyright 2010 The MathWorks, Inc.

% Define the default outputs.
success   = true;
exception = MException.empty;

indx = 1;
specObjs = this.specobjs;

% Loop over all the specification objects and ask them to validate.
while success && indx <= numel(specObjs)
    [success, exception] = validate(specObjs(indx), d);
    indx = indx+1;
end

if nargout
    varargout = {success, exception};
elseif ~success
    throw(exception);
end

% [EOF]
