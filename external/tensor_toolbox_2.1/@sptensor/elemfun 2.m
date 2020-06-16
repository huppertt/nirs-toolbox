function a = elemfun(a,fun)
%ELEMFUN Manipulate the nonzero elements of a sparse tensor.
%
%   X = ELEMFUN(X,@FUN) modifies the elements of X according to the
%   function @FUN which should take and array and output an equally
%   sized array.
%
%   Examples
%   X = sptenrand([10,10,10],10);
%   X = elemfun(X,@sqrt) %<-- square root of every entry
%   X = elemfun(X, @(x) x+1) %<-- increase every entry by 1
%   X = elemfun(X, @(x) x ~= 0) %<-- change every nonzero to be 1
%
%   See also SPTENSOR, SPFUN.
%
%Copyright 2006, Tamara Kolda and Brett Bader, Sandia National Labs


a.vals = fun(a.vals);