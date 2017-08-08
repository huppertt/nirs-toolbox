function help(this, designmethod)
%HELP   Provide help for the specified design method.

%   Author(s): J. Schickler
%   Copyright 2005 The MathWorks, Inc.

hfmethod = feval(getdesignobj(this, designmethod));

help(hfmethod);

% [EOF]
