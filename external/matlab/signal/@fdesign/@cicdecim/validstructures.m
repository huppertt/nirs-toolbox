function vs = validstructures(~, ~, varargin)
%VALIDSTRUCTURES   

%   Copyright 2005-2011 The MathWorks, Inc.

vs_str = {'cicdecim'};

if nargin < 2
    vs.design = vs_str;
else
    vs = vs_str;
end

% [EOF]
