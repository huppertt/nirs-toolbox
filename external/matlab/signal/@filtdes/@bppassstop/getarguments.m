function varargout = getarguments(h, d)
%GETARGUMENTS Returns the standard input arguments

%   Author(s): J. Schickler
%   Copyright 1988-2002 The MathWorks, Inc.

% Get frequency specs, they have been prenormalized
Fstop1 = get(d,'Fstop1');
Fpass1 = get(d,'Fpass1');
Fpass2 = get(d,'Fpass2');
Fstop2 = get(d,'Fstop2');

% Get weights
Wstop1 = get(d,'Wstop1');
Wpass  = get(d,'Wpass');
Wstop2 = get(d,'Wstop2');

args = {[0 Fstop1 Fpass1 Fpass2 Fstop2 1], [0 0 1 1 0 0], [Wstop1 Wpass Wstop2]};

if nargout == 1,
    varargout = {args};
else
    varargout = {args{:}, {}};
end

% [EOF]
