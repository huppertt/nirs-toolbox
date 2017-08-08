function win = generatewindow(this)
%GENERATEWINDOW Generate the window used for the design.

%   Author(s): R. Losada, J. Schickler
%   Copyright 1988-2003 The MathWorks, Inc.

N = get(this,'order');

winobj = get(this, 'WindowObject');

set(winobj, 'Length', N+1);
win = generate(winobj);

% [EOF]
