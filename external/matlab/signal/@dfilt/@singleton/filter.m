function y = filter(Hd,x,dim)
%FILTER Discrete-time filter.
%   Y = FILTER(Hd,X) filters the data X using the discrete-time filter
%   object Hd to create the filtered data Y. The final conditions are
%   stored in Hd.States. 
%  
%   If Hd.PersistentMemory is false (default), initial conditions
%   are set to zero before filtering. 
%  
%   To use non-zero initial conditions, set Hd.PersistentMemory to true
%   and set Hd.States  with a vector of NSTATES(Hd) elements. If a scalar
%   is specified, it will be expanded to a vector of the correct length. 
%
%   FILTER(Hd,X,DIM) operates along the dimension DIM. If X is a vector or 
%   matrix and DIM is 1, every column of X is treated as a channel. If DIM
%   is 2, every row represents a channel.
%
%   See also DFILT.

%   Author: Thomas A. Bryan
%   Copyright 1988-2005 The MathWorks, Inc.

error(nargchk(1,3,nargin,'struct'));

if nargin<2, x = []; end
if nargin<3, dim = []; end

if isempty(x), 
  y = x;
  return; 
end

% Because this filter method used to allow for specifying the states as the
% third input, we check that the DIM input is not a vector/matrix.
if any(size(dim)>1),
    error(message('signal:dfilt:singleton:filter:DimMustBeInt'));
end

% Call super's method
y = super_filter(Hd,x,dim);

