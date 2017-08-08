function [params, values, descs, iargs] = genmcode(h, d)
%GENMCODE Generate MATLAB code

%   Author(s): J. Schickler
%   Copyright 1988-2003 The MathWorks, Inc.

[params, values, descs]  = abstract_genmcode(h,d);

iargs = sprintf('F%s, A, W', getfsstr(d));

% [EOF]
