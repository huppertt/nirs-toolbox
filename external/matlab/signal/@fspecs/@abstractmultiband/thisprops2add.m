function p = thisprops2add(this,varargin)
%THISPROPS2ADD   

%   Author(s): V. Pellissier
%   Copyright 2005 The MathWorks, Inc.

p = propstoadd(this);
% Remove the NormalizedFrequency and Fs properties.
p(1:2) = [];

idx = strmatch('NBands',p);
if ~isempty(idx) && length(varargin)>1,
    % Set NBands first and update propstoadd
    set(this, p{idx}, varargin{idx});
    p = propstoadd(this);
    p(1:2) = [];
end


% [EOF]
