function data = convert2db(this, data)
%CONVERT2DB   COnvert the data to dB.

%   Author(s): J. Schickler
%   Copyright 1988-2004 The MathWorks, Inc.

if nargin < 2
    data = get(this, 'Data');
end

ws = warning; % Cache warning state
warning off   % Avoid "Log of zero" warnings
data = db(abs(data), 'power');  % Call the Convert to decibels engine
warning(ws);  % Reset warning state

% [EOF]
