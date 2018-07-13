function B = repelem(A, varargin)
%REPELEM Replicate elements of a table.
%   B = REPELEM(A, M, N), with table A and scalars M and N, replicates
%   each element of A into an M-by-N block of table B.
%
%   B = REPELEM(A, M, N), with table A and vectors M and N, replicates
%   element A(i, j) into an M(i)-by-N(j) block of table B.
%
%   See also REPMAT, BSXFUN, MESHGRID.

%   Copyright 2014 The MathWorks, Inc.

if nargin ~= 3
    error(message('MATLAB:table:repelem:WrongRHS'));
else
    B = matlab.internal.builtinhelper.repelem(A, varargin{1}, varargin{2});
end
