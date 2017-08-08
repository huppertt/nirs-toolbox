function f = thisisfir(this)
%THISISFIR   Dispatch and call the method.

%   Author(s): J. Schickler
%   Copyright 1988-2004 The MathWorks, Inc.

[b,a] = tf(this);  
f = signalpolyutils('isfir',b,a);

% [EOF]
