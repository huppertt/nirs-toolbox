function varargout = dfgetdistributions(varargin)
%DFGETDISTRIBUTIONS A helper function for dffig2m to get structure defining
%the distributions supported by dfittool. 
%   This internal function is used to avoid generating code that that
%   calls dfswitchyard.

%   Copyright 2012 The MathWorks, Inc.


[varargout{1:nargout}] = dfswitchyard('dfgetdistributions',varargin{:});
