function [params, values, descs, iargs] = genmcode(h, d)
%GENMCODE Generate MATLAB code that designs a lowpass filter using GREMEZ

%   Author(s): J. Schickler
%   Copyright 1988-2003 The MathWorks, Inc.

[params, values, descs]  = abstract_genmcode(h, d);

[fs, fsstr] = getfsstr(d);

iargs = sprintf('[0 Fpass1 Fstop1 Fstop2 Fpass2 %s]%s, %s, %s', fsstr, fs, ...
    '[1 1 0 0 1 1]', '[Dpass1 Dstop Dpass2]');

% [EOF]
