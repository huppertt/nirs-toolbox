function m = hiddenmethods(this)
%HIDDENMETHODS   Return the hidden methods.

%   Author(s): J. Schickler
%   Copyright 2006 The MathWorks, Inc.

% Return the hidden methods of the current specs object.  These methods are
% hidden from the DESIGNMETHODS method, but are accessible if you know
% their names for backwards compatibility reasons.
m = hiddendesigns(this.CurrentSpecs);

% [EOF]
