function b = islphpreorder(this)
%ISLPHPREORDER True filter response is lowpass or highpass

%   Copyright 2008 The MathWorks, Inc.

b = this.F0==0 | this.F0==1;

% [EOF]
