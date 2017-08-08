function [varargout]=insertnan(wasnan,varargin)
%INSERTNAN Insert NaN, space, '' or undefined value into inputs.
%   X1 = internal.stats.insertnan(WASNAN, Y1) inserts missing values in Y1
%   and returns it as X1. WASNAN is a logical column vector and the output
%   of internal.stats.removenan. Its TRUE values indicate the rows in X1
%   that will contain missing values. Y1 is a column vector or a matrix.
%   The type of Y1 can be:
%   Categorical       - X1 is categorical, undefined values represents
%                       missing values.
%   Double            - X1 is double. NaN values represents missing values.
%   Single            - X1 is single. NaN values represents missing values.
%   Character matrix  - X1 is a character matrix. Space represents missing
%                       values.
%   Cell              - X1 is a cell array. empty string '' represents
%                       missing values.
%
%  [X1,X2,...] = internal.stats.insertnan(WASNAN,Y1,Y2,...) accepts any
%  number of input variables Y1,Y2,Y3,.... STATINSERTNAN inserts missing
%  values in Y1, Y2,...  and returns them as X1, X2,... respectively.
%
%  See also REMOVENAN.

%   Copyright 1993-2014 The MathWorks, Inc.

if ~any(wasnan)
     varargout = varargin;
     return;
end

ok = ~wasnan;
len = length(wasnan);
varargout = cell(1,nargin-1);
for j=1:nargin-1
    y = varargin{j};
    if (size(y,1)==1) && sum(ok) > 1
        y =  y';
    end
    
    [~,p] = size(y);
    
    if ischar(y)
        x = repmat(' ', [len,p]);
    elseif isa(y,'categorical')
        x = y([]); % preserve the subclass and attributes
        x(1:len,1:p) = categorical.undefLabel;
    elseif iscell(y)
        x = repmat({''},[len,p]);
    elseif isfloat(y)
        x = nan([len,p],class(y));
    elseif islogical(y)
        error(message('stats:statinsertnan:LogicalInput'));
    else
        error(message('stats:statinsertnan:InputTypeIncorrect'));
    end
    
    x(ok,:) = y;
    
    varargout{j} = x;
end
