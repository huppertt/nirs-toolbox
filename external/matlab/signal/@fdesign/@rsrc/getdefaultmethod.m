function defaultmethod = getdefaultmethod(this)
%GETDEFAULTMETHOD   Get the defaultmethod.

%   Author(s): J. Schickler
%   Copyright 2005 The MathWorks, Inc.

if strcmpi(this.Response,'nyquist'),
    defaultmethod = 'kaiserwin';
else
    defaultmethod = 'equiripple';
end

% [EOF]
