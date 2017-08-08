function b = isdesignmethod(this, method)
%ISDESIGNMETHOD   Returns true if the method is a valid designmethod.

%   Author(s): J. Schickler
%   Copyright 1999-2005 The MathWorks, Inc.

d = designmethods(this);
if isa(method, 'function_handle'),
    method = func2str(method);
end

b = any(strcmpi(method, d));

% If the method is not listed it might be a "hidden" design.  These are
% usually designs that we no longer document and do not want to list, but
% we still want ISDESIGNMETHODS to return true for backwards compatibility.
if ~b
    hdesigns = hiddendesigns(this.CurrentSpecs);
    b = any(strcmpi(method, hdesigns));
end

% [EOF]
