function Hd = dispatch(this)
%DISPATCH   Return the lwdfilt.  

%   Author(s): V. Pellissier
%   Copyright 2005 The MathWorks, Inc.

Hd = lwdfilt.tf([zeros(1,this.Latency) 1]);

% [EOF]
