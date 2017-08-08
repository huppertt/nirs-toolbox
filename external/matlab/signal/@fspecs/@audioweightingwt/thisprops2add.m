function p = thisprops2add(this,varargin)
%THISPROPS2ADD   

%   Copyright 2009 The MathWorks, Inc.

p = fieldnames(this);

% Remove the ResponseType, NormalizedFrequency and Fs properties.
p(1:3) = [];

% Remove WeightingType
p(end) = [];

% [EOF]
