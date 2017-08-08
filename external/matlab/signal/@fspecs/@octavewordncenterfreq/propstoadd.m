function p = propstoadd(this,varargin)
%PROPSTOADD   

%   Author(s): V. Pellissier
%   Copyright 2005 The MathWorks, Inc.

p = fieldnames(this);

% Remove the ResponseType
p(1) = [];

% Remove privLthOctave
p(end) = [];

% [EOF]
