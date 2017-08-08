function specs = whichspecs(h)
%WHICHSPECS Determine which specs are required for this class.

%   Author(s): J. Schickler
%   Copyright 1988-2002 The MathWorks, Inc.

% Call super's method
specs = fps_whichspecs(h);

specs(1).defval = 2880;
specs(2).defval = 3360;

% [EOF]
