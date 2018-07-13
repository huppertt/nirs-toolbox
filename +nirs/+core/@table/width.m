function w = width(t)
%WIDTH Number of variables in a table.
%   W = WIDTH(T) returns the number of variables in the table T.  WIDTH(T) is
%   equivalent to SIZE(T,2).  Note that variables in a table may themselves have
%   multiple columns.  WIDTH(T) does not account for that.
%  
%   See also HEIGHT, SIZE, NUMEL.

%   Copyright 2012 The MathWorks, Inc. 

w = t.nvars;
