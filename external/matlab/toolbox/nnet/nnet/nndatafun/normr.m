function y = normr(x)
%NORMR Normalize rows of matrices.
%
%  <a href="matlab:doc normr">normr</a>(X) takes a single matrix or cell array of matrices and returns
%  the matrices with rows normalized to a length of one.
%
%  Here the rows of a random matrix are normalized.
%
%    x = <a href="matlab:doc rands">rands</a>(4,8);
%    y = <a href="matlab:doc normr">normr</a>(x)
%
%  See also NORMC.

% Copyright 1992-2015 The MathWorks, Inc.

% Checks
if nargin < 1
    error(message('nnet:Args:NotEnough')); 
end
wasMatrix = ~iscell(x);
x = nntype.data('format',x,'Data');

% Compute
y = cell(size(x));
for i=1:numel(x)
    xi = x{i};
    xi(~isfinite(xi)) = 0;
    len = sqrt(sum(xi.^2,2));
    yi = bsxfun(@rdivide,xi,len);
    zeroRows = find(len==0);
    if ~isempty(zeroRows)
        numColumns = size(xi,2);
        row = ones(1,numColumns) ./ sqrt(numColumns);
        yi(zeroRows,:) = repmat(row,numel(zeroRows),1);
    end
    y{i} = yi;
end

% Format
if wasMatrix
    y = y{1};
end
