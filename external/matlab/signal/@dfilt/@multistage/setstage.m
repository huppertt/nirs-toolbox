function setstage(this, Hd, pos)
%SETSTAGE   Set the stage.

%   Author(s): J. Schickler
%   Copyright 1988-2004 The MathWorks, Inc.

error(nargchk(3,3,nargin,'struct'));

s = this.Stage;
s(pos) = Hd;
this.Stage = s;

% [EOF]
