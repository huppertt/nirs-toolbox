function varargout = getarguments(h, d)
%GETARGUMENTS Returns the standard input arguments

%   Author(s): J. Schickler
%   Copyright 1988-2002 The MathWorks, Inc.

% Get frequency specs, they have been prenormalized
Fpass1 = get(d,'Fpass1');
Fstop1 = get(d,'Fstop1');
Fstop2 = get(d,'Fstop2');
Fpass2 = get(d,'Fpass2');

% Get weights
Wpass1 = get(d,'Wpass1');
Wstop  = get(d,'Wstop');
Wpass2 = get(d,'Wpass2');

args = {[0 Fpass1 Fstop1 Fstop2 Fpass2 1], [1 1 0 0 1 1], [Wpass1 Wstop Wpass2]};

if nargout == 1,
    varargout = {args};
else
    varargout = {args{:}, {}};
end

% [EOF]
