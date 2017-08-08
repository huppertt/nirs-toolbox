function [N,F,E,A,nfpts,Fs,normFreqFlag] = validatespecs(this)
%VALIDATESPECS Validate the specs

%   Copyright 2011 The MathWorks, Inc.

% Get filter order, amplitudes and frequencies 
N = this.FilterOrder;
[F,E,A,nfpts,Fs,normFreqFlag] = super_validatespecs(this);

% [EOF]
