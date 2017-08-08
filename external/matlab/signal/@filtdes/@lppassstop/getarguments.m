function varargout = getarguments(h, d)
%GETARGUMENTS Returns the standard input arguments

%   Author(s): J. Schickler
%   Copyright 1988-2002 The MathWorks, Inc.

% Get frequency specs, they have been prenormalized
Fpass = get(d,'Fpass');
Fstop = get(d,'Fstop');

% Get weights
Wpass = get(d,'Wpass');
Wstop = get(d,'Wstop');

args = {[0 Fpass Fstop 1], [1 1 0 0], [Wpass Wstop]};

if nargout == 1,
    varargout = {args};
else
    varargout = {args{:}, {}};
end

% [EOF]
