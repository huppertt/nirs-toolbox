function [params, values, descs, iargs] = genmcode(h, d)
%GENMCODE Generate MATLAB code that designs a lowpass filter using GREMEZ

%   Author(s): J. Schickler
%   Copyright 1988-2003 The MathWorks, Inc.

[params, values, descs] = abstract_genmcode(h, d);

[fs, fsstr] = getfsstr(d);

iargs = sprintf('[0 Fpass Fstop %s]%s, [1 1 0 0], [Dpass Dstop]', fsstr, fs);

% [EOF]
