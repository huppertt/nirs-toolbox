function disp(this)
%DISP   Object Display.

%   Author(s): J. Schickler
%   Copyright 1988-2004 The MathWorks, Inc.

s = dispstr(this);

f1 = min(find(s == char(10)));

s = [blanks(4) s(1:f1) blanks(4) s(f1+1:end) char(10)];

disp(s);

% [EOF]
