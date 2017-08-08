function cmd = maskinfo(hObj, d)
%MASKINFO Return the mask information

%   Author(s): J. Schickler
%   Copyright 1988-2003 The MathWorks, Inc.

cmd{1}.frequency  = get(d, 'FrequencyVector');
cmd{1}.freqfcn    = 'aline';
cmd{1}.amplitude  = ones(length(d.FrequencyVector), 1);
cmd{1}.properties = {'Color', 'red', 'linestyle', '--'};

% [EOF]
