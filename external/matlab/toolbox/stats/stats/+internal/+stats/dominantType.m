function proto = dominantType(varargin)
%DOMINANTTYPE return a prototype of the dominant input type
%   DOMINANTTYPE uses standard arithmetic rules to determine a prototype of
%   the dominant type from the inputs.
%
%   Example:
%   dominantType( single(1), double(2) ) => single
%   dominantType( single(1), gpuArray(2) ) => gpuArray(single)

%   Copyright 2014 The MathWorks, Inc.

% Add together empty arrays so that the types propagate but no actual
% calculation is required.
proto = varargin{1}([]);
for ii=2:nargin
    proto = proto + varargin{ii}([]);
end
