function [x,perm,nshifts] = shiftdata(x,dim)
%SHIFTDATA  Shift data to operate on a specified dimension.
%   [X,PERM,NSHIFTS] = SHIFTDATA(X,DIM) shifts data X so that dimension DIM is
%   permuted to the first column using the same permutation as the built-in
%   FILTER function.  The permutation vector that is used is returned in
%   vector PERM.
%
%   If DIM is missing or empty, then the first non-singleton dimension is
%   shifted to the first column, and the number of shifts is returned in
%   NSHIFTS. 
%  
%   SHIFTDATA is handy for creating functions that work along a certain
%   dimension, like FILTER, GOERTZEL, SGOLAYFILT and SOSFILT.
%
%   SHIFTDATA is meant to be used in tandem with UNSHIFTDATA, which shifts
%   the data back to the original shape.
%
%   Example:
%     x = magic(3)
%     [x,perm,nshifts] = shiftdata(x,2) % Work along 2nd dimension
%     y = unshiftdata(x,perm,nshifts)   % Reshapes back to original
%
%     x = 1:5                            % Originally a row
%     [x,perm,nshifts] = shiftdata(x,[]) % Work along 1st nonsingleton dimension
%     y = unshiftdata(x,perm,nshifts)    % Reshapes back to original
%
%   See also UNSHIFTDATA, PERMUTE, SHIFTDIM.
  
%   Copyright 1988-2013 The MathWorks, Inc.

if nargin<2
  dim=[];
end

if isempty(dim)
  % Work along the first nonsingleton dimension
  [x,nshifts] = shiftdim(x);
  perm = [];
else
  % Put DIM in the first dimension (this matches the order that the built-in
  % filter uses in function attfcn:mfDoDimsPerm)
  
  % Cast to enforce precision rules
  dim = signal.internal.sigcasttofloat(dim,'double','shiftdata','DIM',...
    'allownumeric');
   
  perm = [dim,1:dim-1,dim+1:ndims(x)];
  x = permute(x,perm);
  nshifts = [];
end
