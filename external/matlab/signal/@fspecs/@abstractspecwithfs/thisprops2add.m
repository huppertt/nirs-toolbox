function p = thisprops2add(this,varargin)
%THISPROPS2ADD   

%   Author(s): V. Pellissier
%   Copyright 2005 The MathWorks, Inc.

p = propstoadd(this);

% Remove the NormalizedFrequency and Fs properties.
p(1:2) = [];


% [EOF]
