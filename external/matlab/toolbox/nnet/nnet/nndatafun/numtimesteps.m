function ts = numtimesteps(x)
%NUMTIMESTEPS Number of samples in neural network data.
%
%  <a href="matlab:doc numtimesteps">numtimesteps</a>(X) returns the number of timesteps in X, which must be
%  NN data in either matrix or cell array form.
%
%  If X is a matrix the result is 1.
%
%  If X is a cell array the result is the number of rows in X.
%
%  This code calculates the number of timesteps represented by matrix data:
%
%    x = [1 2 3; 4 7 4]
%    n = <a href="matlab:doc numtimesteps">numtimesteps</a>(x)
%
%  This code calculates the number of timesteps represented by cell data:
%
%    x = {[1:3; 4:6] [7:9; 10:12]; [13:15] [16:18]}
%    n = numtimesteps(x)
%
%  See also GETTIMESTEPS, SETTIMESTEPS, CATTIMESTEPS, NNDATA, NNSIZE.

% Copyright 2010 The MathWorks, Inc.

if nargin < 1,error(message('nnet:Args:NotEnough')); end
x = nntype.data('format',x,'Data');

if iscell(x)
  ts = nnfast.numtimesteps(x);
else
  ts = 1;
end
