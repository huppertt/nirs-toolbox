function HdB = convert2db(this,H)
%CONVERT2DB   Convert input response to db values.

%   Author(s): P. Pacheco
%   Copyright 1988-2003 The MathWorks, Inc.

ws = warning; % Cache warning state
warning off   % Avoid "Log of zero" warnings
HdB = db(H,'voltage');  % Call the Convert to decibels engine
warning(ws);  % Reset warning state

% [EOF]
