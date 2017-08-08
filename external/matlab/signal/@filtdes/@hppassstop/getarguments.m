function varargout = getarguments(h, d)
%GETARGUMENTS Returns the standard input arguments

%   Author(s): J. Schickler
%   Copyright 1988-2002 The MathWorks, Inc.

% Get frequency specs, they have been prenormalized
Fstop = get(d,'Fstop');
Fpass = get(d,'Fpass');

% Get weights
Wstop = get(d,'Wstop');
Wpass = get(d,'Wpass');

args = {[0 Fstop Fpass 1], [0 0 1 1], [Wstop Wpass]};

if nargout == 1,
    varargout = {args};
else
    varargout = {args{:}, {}};
end    

% [EOF]
