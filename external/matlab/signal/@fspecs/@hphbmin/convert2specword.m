function [fspecsword,dm,dopts,Nstep] = convert2specword(this,cfmethod,N)
%CONVERT2SPECWORD Convert minimum order spec to spec with order for
%                 constrained worrd length FIR design

% This should be a private method

%   Copyright 2008 The MathWorks, Inc.

if strcmpi(cfmethod.DesignAlgorithm,'Equiripple'),
    dm = 'equiripple';
    dopts = {'FilterStructure', cfmethod.FilterStructure, ...
        'MinPhase', cfmethod.MinPhase, ...
        'StopbandShape', cfmethod.StopbandShape, ...
        'StopbandDecay', cfmethod.StopbandDecay};
    
elseif strcmpi(cfmethod.DesignAlgorithm,'Kaiser window'),
    dm = 'kaiserwin';
    dopts = {'FilterStructure', cfmethod.FilterStructure};
else
    error(message('signal:fspecs:hphbmin:convert2specword:InvalidMethod'));
end

fspecsword = fspecs.hphbordntw;
fspecsword.FilterOrder = N;
if ~this.NormalizedFrequency,
    normalizefreq(fspecsword,false,this.Fs)
end
fspecsword.TransitionWidth = this.TransitionWidth;

% Filter Order increment
Nstep = 4;

% [EOF]
