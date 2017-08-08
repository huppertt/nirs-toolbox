function [p, v, d, ia] = ifir_genmcode(h)
%IFIR_GENMCODE Return the IFIR code

%   Author(s): J. Schickler
%   Copyright 1988-2002 The MathWorks, Inc.

p = {'L'};
v = {getmcode(h, 'InterpolationFactor')};
d = {'Interpolation Factor'};

opt = get(h, 'Optimization');

if strcmpi(opt, 'intermediate'),
    ia = '';
else
    ia = ', optim';
    p  = {p{:}, 'optim'};
    v  = {v{:}, sprintf('''%s''', opt)};
    d  = {d{:}, 'Optimization Level'};
end

% [EOF]
