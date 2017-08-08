function hresp = responseobj(this)
%RESPONSEOBJ   Magresp response object.
%
% This is a private method.

%   Author(s): P. Pacheco
%   Copyright 1988-2006 The MathWorks, Inc.

% Create the response object. 
hresp = sigresp.freqz(this);

% [EOF]
