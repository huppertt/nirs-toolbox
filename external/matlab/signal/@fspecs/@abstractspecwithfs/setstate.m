function setstate(this, state)
%SETSTATE   Set the state of the object.

%   Author(s): J. Schickler
%   Copyright 1988-2003 The MathWorks, Inc.

set(this, rmfield(state, 'ResponseType'));

% [EOF]
