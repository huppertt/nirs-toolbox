function [badin,wasnan,varargout]=removenan(varargin)
%REMOVENAN Remove NaN values from inputs
%   [~,WASNAN,X1] = internal.stats.removenan(Y1) removes missing values
%   from Y1 and returns it as X1. WASNAN is a logical column vector
%   indicating where there were NaN values removed. Y1 is a column vector
%   or a matrix. The type of Y1 can be:
%
%   Categorical       - X1 is categorical, undefined values represents
%                       missing values.
%   Double            - X1 is double. NaN values represents missing values.
%   Single            - X1 is single. NaN values represents missing values.
%   Character matrix  - X1 is a character matrix. Space represents missing
%                       values.
%   Cell              - X1 is a cell array. empty string '' represents
%                       missing values.
%
%  [BAD,WASNAN,X1,x2,...] = internal.stats.removenan(Y1,Y2,...) accepts any
%  number of input variables Y1,Y2,Y3,..., removes missing values, and
%  returns them as X1, X2,... respectively. All variables should have the
%  same number of observations, and BAD is 0 if that is true. If they to
%  not have the same number, BAD returns the index of a variable whose
%  length does not match the others.
%
%  See also INSERTNAN.

%   Copyright 1993-2014 The MathWorks, Inc.

varargoutsize = nargout - 2;
if varargoutsize > 0
    [badin,wasnan,varargout{1:varargoutsize}] = statslib.internal.removenan(varargin{:});
else
    [badin,wasnan] = statslib.internal.removenan(varargin{:});
end
