function [p,v,d,args] = cremez_genmcode(h)
%CREMEZ_GENMCODE Get the CREMEZ MATLAB code information

%   Author(s): J. Schickler
%   Copyright 1988-2002 The MathWorks, Inc.

sym  = get(h, 'SymmetryConstraint');
if ~strcmpi(sym, 'default')
    args = sprintf(', symc');
    p = {'symc'};
    v = {sprintf('''%s''', lower(sym))};
    d = {'Symmetry Constraint'};
else
    args = '';
    p = {};
    v = {};
    d = {};
end

optim = get(h, 'SecondStageOptimization');
if strcmpi(optim, 'off')
    args = sprintf('%s, ''skip_stage2''', args);
end

p = {p{:}, 'debug', 'dens'};
v = {v{:}, sprintf('''%s''', lower(get(h, 'DebugMode'))), ...
        sprintf('%d', get(h, 'DensityFactor'))};
d = {d{:}, 'Debug Mode', ''};

args = sprintf('%s, {dens}, debug', args);

% [EOF]
