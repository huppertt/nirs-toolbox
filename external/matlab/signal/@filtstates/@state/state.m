function this = state(value)
%STATE   Construct a STATE object.

%   Author(s): J. Schickler
%   Copyright 1988-2004 The MathWorks, Inc.

this = filtstates.state;

if nargin
    this.Value = value;
end

% [EOF]
