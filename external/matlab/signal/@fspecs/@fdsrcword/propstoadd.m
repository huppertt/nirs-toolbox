function p = propstoadd(this)
%PROPSTOADD   

%   Copyright 2007 The MathWorks, Inc.

p = fieldnames(this);

% Remove the ResponseType
p(1) = [];

% Remove privInterpolationFactor and privDecimationFactor
p(end-1:end) = [];

