function chknostropts(varargin)
%CHKNOSTROPTS Checks that no strings are specified
%   errors on the first encountered string
%
%   This file is for internal use only

%   Copyright 2013 The MathWorks, Inc.

strOpts = varargin(cellfun(@ischar, varargin));
if ~isempty(strOpts)
  ME = MException(message('signal:chknostropts:UnknownOption',strOpts{1}));
  throwAsCaller(ME);
end

