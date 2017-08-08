function [params, values, descs, iargs] = genmcode(h, d)
%GENMCODE Generate MATLAB code that designs a lowpass filter using GREMEZ

%   Author(s): J. Schickler
%   Copyright 1988-2003 The MathWorks, Inc.

mu = get(d, 'magUnits'); set(d, 'magUnits', 'linear');

params = {'Fpass1', 'Fstop1', 'Fstop2', 'Fpass2', 'Wpass1', 'Wstop', 'Dpass2'};
values = {getmcode(d, 'Fpass1'), getmcode(d, 'Fstop1'), getmcode(d, 'Fstop2'), ...
        getmcode(d, 'Fpass2'), getmcode(d, 'Wpass1'), getmcode(d, 'Wstop'), ...
        getmcode(d, 'Dpass2')};
descs  = {'', '', '', '', '', '', ''};

set(d, 'magUnits', mu);

[fsstr, fs] = getfsstr(d);

iargs = sprintf('[0 Fpass1 Fstop1 Fstop2 Fpass2 %s]%s, %s, %s', ...
    fs, fsstr, '[1 1 0 0 1 1]', '[Wpass1 Wstop Dpass2]');

% [EOF]
