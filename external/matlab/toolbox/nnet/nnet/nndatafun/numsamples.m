function q = numsamples(x)
%NUMSAMPLES Number of samples in neural network data.
%
%  <a href="matlab:doc numsamples">numsamples</a>(X) return the number of samples of X, which must be NN data
%  in either matrix or cell array form.
%
%  If X is a matrix the result is the number of columns.
%
%  If X is a cell array the result is the number of columns of the
%  matrices in X.
%
%  This code calculates the number of samples represented by matrix data:
%
%    x = [1 2 3; 4 7 4]
%    n = <a href="matlab:doc numsamples">numsamples</a>(x)
%
%  This code calculates the number of samples represented by cell data:
%
%    x = {[1:3; 4:6] [7:9; 10:12]; [13:15] [16:18]}
%    n = <a href="matlab:doc numsamples">numsamples</a>(x)
%
%  See also SETSAMPLES, GETSAMPLES, CATSAMPLES, ISNNDATA, NNDIMENSIONS

% Copyright 2010 The MathWorks, Inc.

if nargin < 1,error(message('nnet:Args:NotEnough')); end
x = nntype.data('format',x,'X');

if iscell(x)
  q = nnfast.numsamples(x);
else
  q = size(x,2);
end
