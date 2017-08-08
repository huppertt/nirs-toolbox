function [N,F,E,A,nfpts,Fs,NormFreqFlag] = validatespecs(this)
%VALIDATESPECS Validate the specs

%   Copyright 2005-2011 The MathWorks, Inc.

% Get filter order, amplitudes and frequencies 
N = this.FilterOrder;
[F,E,A,nfpts,Fs,NormFreqFlag] = super_validatespecs(this);

% [EOF]
