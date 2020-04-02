function t = cat(dim,varargin)
%CAT Concatenate tables.
%   T = CAT(DIM, T1, T2, ...) concatenates the tables T1, T2, ... along
%   dimension DIM by calling the TABLE/HORZCAT or TABLE/VERTCAT method.
%   DIM must be 1 or 2.
%
%   See also HORZCAT, VERTCAT.

%   Copyright 2012 The MathWorks, Inc.

if dim == 1
    t = vertcat(varargin{:});
elseif dim == 2
    t = horzcat(varargin{:});
else
    error(message('MATLAB:table:cat:InvalidDim'));
end
