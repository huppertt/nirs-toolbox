function p = thisprops2add(this,varargin)
%THISPROPS2ADD   

%   Author(s): V. Pellissier
%   Copyright 2005-2006 The MathWorks, Inc.

p = fieldnames(this);

% Remove the ResponseType, NormalizedFrequency and Fs properties.
p(1:3) = [];

% Remove privFracdelay
p(end) = [];

% [EOF]
