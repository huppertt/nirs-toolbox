function [params, values, descs, str, args] = genmcode(h, d)
%GENMCODE Generate MATLAB code

%   Author(s): J. Schickler
%   Copyright 1988-2002 The MathWorks, Inc.

[fsstr, fs] = getfsstr(d);

[params, values, descs] = abstract_genmcode(h, d);

str = sprintf('\nF = [0 Fstop Fpass %s]%s;\n', fs, fsstr);

args = 'F, F, [0 0 1 1], [Wstop Wstop Wpass Wpass]';

% [EOF]
