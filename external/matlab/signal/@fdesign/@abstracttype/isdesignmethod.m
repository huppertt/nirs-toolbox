function b = isdesignmethod(this, method)
%ISDESIGNMETHOD   Returns true if the method is a valid designmethod.

%   Author(s): J. Schickler
%   Copyright 1999-2004 The MathWorks, Inc.

d = designmethods(this);

if isa(method, 'function_handle'),
    method = func2str(method);
end

b = any(strcmpi(method, d));

% [EOF]
