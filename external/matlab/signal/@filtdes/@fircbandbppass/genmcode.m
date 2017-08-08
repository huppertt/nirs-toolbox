function [params, values, descs, iargs] = genmcode(h, d)
%GENMCODE Generate MATLAB code that designs a lowpass filter using GREMEZ

%   Author(s): J. Schickler
%   Copyright 1988-2003 The MathWorks, Inc.

mu = get(d, 'magUnits'); set(d, 'magUnits', 'linear');

params = {'N', 'Fstop1', 'Fpass1', 'Fpass2', 'Fstop2', 'Wstop1', 'Dpass', 'Wstop2'};
values = {getmcode(d, 'Order'), getmcode(d, 'Fstop1'), getmcode(d, 'Fpass1'), ...
        getmcode(d, 'Fpass2'), getmcode(d, 'Fstop2'), getmcode(d, 'Wstop1'), ...
        getmcode(d, 'Dpass'), getmcode(d, 'Wstop2')};
descs  = {'', '', '', '', '', '', '', ''};

set(d, 'magUnits', mu);

[fsstr, fs] = getfsstr(d);

iargs = sprintf('[0 Fstop1 Fpass1 Fpass2 Fstop2 %s]%s, %s, %s]', ...
    fs, fsstr, '[0 0 1 1 0 0]', '[Wstop1 Dpass Wstop2]');

% [EOF]
