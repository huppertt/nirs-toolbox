function [params, values, descs, str, args] = genmcode(h, d)
%GENMCODE Generate MATLAB code

%   Author(s): J. Schickler
%   Copyright 1988-2002 The MathWorks, Inc.

fsstr = getfsstr(d);

params = {'F', 'E', 'A', 'W'};
values = {getmcode(d,'FrequencyVector'), getmcode(d,'FrequencyEdges'), ...
        getmcode(d,'MagnitudeVector'), getmcode(d,'WeightVector')};
descs  = {'', '', '', ''};

str = '';

args = sprintf('F%s, E%s, A, W', fsstr, fsstr);

% [EOF]
