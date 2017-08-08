function [varargout]=statinsertnan(wasnan,varargin)
%STATINSERTNAN Insert NaN, space, '' or undefined value into inputs.
%   X1 = STATINSERTNAN(WASNAN, Y1) inserts missing values in Y1 and returns
%   it as X1. WASNAN is a logical column vector and the output of
%   STATREMOVENAN. Its TRUE values indicate the rows in X1 that will
%   contain missing values. Y1 is a column vector or a matrix. The type of
%   Y1 can be:
%   Categorical       - X1 is categorical, undefined values represents
%                       missing values.
%   Double            - X1 is double. NaN values represents missing values.
%   Single            - X1 is single. NaN values represents missing values.
%   Character matrix  - X1 is a character matrix. Space represents missing
%                       values.
%   Cell              - X1 is a cell array. empty string '' represents
%                       missing values.
%
%  [X1,X2,...] = STATINSERTNAN(WASNAN,Y1,Y2,...) accepts any number of
%  input variables Y1,Y2,Y3,.... STATINSERTNAN inserts missing values in
%  Y1, Y2,...  and returns them as X1, X2,... respectively.
%
%  See also STATREMOVENAN.


%   Copyright 1993-2014 The MathWorks, Inc.


[varargout{1:nargout}] = insertnan(wasnan,varargin{:});
