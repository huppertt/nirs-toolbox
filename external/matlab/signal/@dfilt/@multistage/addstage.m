function addstage(this, Hd, pos)
%ADDSTAGE   Add a stage to the filter.

%   Author(s): J. Schickler
%   Copyright 1988-2004 The MathWorks, Inc.

error(nargchk(2,3,nargin,'struct'));

if nargin==2,
    pos = length(this.Stage)+1;
end

this.Stage = [this.Stage(1:pos-1); Hd; this.Stage(pos:end)];

% [EOF]
