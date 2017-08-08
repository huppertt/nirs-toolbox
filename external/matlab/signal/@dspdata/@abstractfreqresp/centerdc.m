function centerdc(this,state)
%CENTERDC   Shift the zero-frequency component to center of spectrum.

%   Author(s): P. Pacheco
%   Copyright 1988-2003 The MathWorks, Inc.

if nargin < 2,
    state = true;
end

if ~xor(this.centerdc, state)
    return;  % State specified = object's state, therefore no-op.
end
this.centerdc = state;

% Call subclasses' method.
thiscenterdc(this);

% [EOF]
