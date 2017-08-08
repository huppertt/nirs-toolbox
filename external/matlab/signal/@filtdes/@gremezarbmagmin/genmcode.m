function [params, values, descs, iargs] = genmcode(h, d)
%GENMCODE Generate MATLAB code

%   Author(s): J. Schickler
%   Copyright 1988-2003 The MathWorks, Inc.

params = {'F', 'A', 'R'};
values = {getmcode(d, 'FrequencyVector'), getmcode(d, 'MagnitudeVector'), ...
        getmcode(d, 'RippleVector')};
descs = {'', '', ''};

iargs = sprintf('F%s, A, R', getfsstr(d));

% [EOF]
