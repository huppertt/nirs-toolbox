function hresp = responseobj(this)
%RESPONSEOBJ   Pseudopowerresp response object.
%
% This is a private method.

%   Author(s): P. Pacheco
%   Copyright 1988-2006 The MathWorks, Inc.

% Create the response obj. 
hresp = sigresp.pseudopowerresp(this);

% [EOF]
