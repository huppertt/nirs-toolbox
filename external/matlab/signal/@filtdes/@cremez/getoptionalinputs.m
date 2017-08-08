function args = getoptionalinputs(d)
%GETOPTIONALINPUTS Get the optional inputs to the design method.

%   Author(s): J. Schickler
%   Copyright 1988-2002 The MathWorks, Inc.

args = {{get(d, 'DensityFactor')}};

sym  = get(d, 'SymmetryConstraint');
if ~strcmpi(sym, 'default')
    args = {args{:}, lower(sym)};
end

optim = get(d, 'SecondStageOptimization');
if strcmpi(optim, 'off')
    args = {args{:}, 'skip_stage2'};
end

debug = get(d, 'DebugMode');
if ~strcmpi(debug, 'off'),
    args = {args{:}, debug};
end

% [EOF]