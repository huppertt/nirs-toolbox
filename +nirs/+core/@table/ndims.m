function n = ndims(t)
%NDIMS Number of dimensions of a table.
%   N = NDIMS(T) returns the number of dimensions in the table T.  The
%   number of dimensions in a table is always 2.
%  
%   See also SIZE.

%   Copyright 2012-2013 The MathWorks, Inc. 

n = t.ndims;
