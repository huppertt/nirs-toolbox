function a = cat(dim,varargin)
%CAT Concatenate dataset arrays.
%   DS = CAT(DIM, DS1, DS2, ...) concatenates the dataset arrays DS1, DS2,
%   ... along dimension DIM by calling the @DATASET/HORZCAT or
%   @DATASET/VERTCAT method. DIM must be 1 or 2.
%
%   See also DATASET/HORZCAT, DATASET/VERTCAT.

%   Copyright 2006-2011 The MathWorks, Inc.


if dim == 1
    a = vertcat(varargin{:});
elseif dim == 2
    a = horzcat(varargin{:});
else
    error(message('stats:dataset:cat:InvalidDim'));
end
