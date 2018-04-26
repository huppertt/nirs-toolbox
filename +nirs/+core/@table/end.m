function e = end(t,k,n)
%END Last index in an indexing expression for a table.
%   END(T,K,N) is called for indexing expressions involving the table T
%   when END is part of the K-th index out of N indices.  For example, the
%   expression T(end-1,:) calls T's END method with END(T,1,2).
%
%   See also SIZE.

%   Copyright 2012 The MathWorks, Inc. 

switch k
case 1
    e = t.nrows;
case 2
    e = t.nvars;
otherwise
    e = 1;
end
