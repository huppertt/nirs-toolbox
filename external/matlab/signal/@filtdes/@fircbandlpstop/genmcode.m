function [params, values, descs, iargs] = genmcode(h, d)
%GENMCODE Generate MATLAB code that designs a lowpass filter using GREMEZ

%   Author(s): J. Schickler
%   Copyright 1988-2003 The MathWorks, Inc.

mu = get(d, 'magUnits'); set(d, 'magUnits', 'linear');

params = {'N', 'Fpass', 'Fstop', 'Dstop'};
values = {getmcode(d, 'Order'), getmcode(d, 'Fpass'), getmcode(d, 'Fstop'), ...
        getmcode(d, 'Dstop')};
descs  = {'', '', '', ''};

set(d, 'magUnits', mu);

[fsstr, fs] = getfsstr(d);

iargs = sprintf('[0 Fpass Fstop %s]%s, [1 1 0 0], [1 Dstop]', fs, fsstr);

% [EOF]
