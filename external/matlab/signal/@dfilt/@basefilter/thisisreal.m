function b = thisisreal(this)
%THISISREAL   Dispatch and call the method.

%   Author(s): J. Schickler
%   Copyright 1988-2010 The MathWorks, Inc.

error(message('signal:dfilt:basefilter:thisisreal:AbstractFunction', class( this )));

% [EOF]
