function [badin,wasnan,varargout]=statremovenan(varargin)
%STATREMOVENAN Remove NaN values from inputs

%   Copyright 1993-2012 The MathWorks, Inc.


[badin,wasnan,varargout{1:nargout-2}] = internal.stats.removenan(varargin{:});
