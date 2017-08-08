function p = propstoadd(this,varargin)
%PROPSTOADD   

%   Copyright 2009 The MathWorks, Inc.

p = fieldnames(this);

% Remove the ResponseType
p(1) = [];

% [EOF]
