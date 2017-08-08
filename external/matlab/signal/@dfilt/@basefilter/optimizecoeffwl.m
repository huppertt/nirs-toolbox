function [Hbest,mrfflag] = optimizecoeffwl(this,varargin)
% This should be a private method.

%   Copyright 2009 The MathWorks, Inc.

% Make sure to work with reference filter in case filter has been quantized
Href = reffilter(this);

try
    [Hbest,mrfflag] = optimizecoeffwlfir(this,Href,varargin{:});
catch ME
    throw(ME);
end




% [EOF]
