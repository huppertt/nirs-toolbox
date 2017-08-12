function n = numelements(x)
%NUMELEMENTS Number of elements in neural network data.
%
%  <a href="matlab:doc numelements">numelements</a>(X) returns the number of elements in each signal of X,
%  which must be NN data in either matrix or cell array form.
%
%  If X is a matrix the result is the number of rows of X.
%
%  If X is a cell array the result is an Sx1 vector, where S is the
%  number of signals (i.e. rows of X), and each element S(i) is the
%  number of elements in each signal i (i.e. rows of X{i,1}.
%
%  This code calculates the number of elements represented by matrix data:
%
%    x = [1 2 3; 4 7 4]
%    n = <a href="matlab:doc numelements">numelements</a>(x)
%
%  This code calculates the number of elements represented by cell data:
%
%    x = {[1:3; 4:6] [7:9; 10:12]; [13:15] [16:18]}
%    n = <a href="matlab:doc numelements">numelements</a>(x)
%
%  See also SETELEMENTS, GETELEMENTS, CATELEMENTS, NNDATA, NNSIZE.

% Copyright 2010 The MathWorks, Inc.

if nargin < 1,error(message('nnet:Args:NotEnough')); end
x = nntype.data('format',x,'Data');

if iscell(x)
  n = nnfast.numelements(x);
else
  n = size(x,1);
end
