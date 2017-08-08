function b = islphpreorder(this)
%ISLPHPREORDER True filter response is lowpass or highpass

%   Copyright 2008 The MathWorks, Inc.

b = this.Flow==0 | this.Fhigh==1;

% [EOF]
