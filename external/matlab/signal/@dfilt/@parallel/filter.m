function y = filter(Hd,x,dim)
%FILTER Discrete-time filter.
%   Y= FILTER(Hd,X,DIM) filters data X over dimension DIM with and returns
%   output Y. 
%
%   See also DFILT.   
  
%   Author: Thomas A. Bryan
%   Copyright 1988-2013 The MathWorks, Inc.

% parallel/filter is not supported for fixed-pt filters
if ~isparallelfilterable(Hd)
    error(message('signal:dfilt:parallel:filter:ParallelFixedFilter'));
end

% Check if all stages have the same overall rate change factor
checkvalidparallel(Hd);

if nargin<3, dim=[]; end

if isempty(x), 
  y = x;
  return; 
end

% Because this filter method used to allow for specifying the states as the
% third input, we check that the DIM input is not a vector/matrix.
if any(size(dim)>1),
    error(message('signal:dfilt:parallel:filter:DimMustBeInt'));
end

% Get current value of reset states
resetval = Hd.PersistentMemory;

y = 0;

flagdiff = false;

needTranspose = any(size(x)==1) && isempty(dim) && isrow(x);
if needTranspose
  x = x(:);
end

for k=1:length(Hd.Stage)
    if Hd.Stage(k).PersistentMemory~=resetval,
        Hd.Stage(k).PersistentMemory = resetval;
        flagdiff = true;
    end
    yk = filter(Hd.Stage(k),x,dim);
    if k==1
      initSize = size(yk);
    elseif ~isequal(size(yk),initSize)
      error(message('signal:dfilt:parallel:filter:MismatchBranches'));
    end
    y = y + yk;
end

if needTranspose
  y = y.';
end

if flagdiff,
     warning(message('signal:dfilt:parallel:filter:flagInconsistency', mat2str( resetval )));
end

% Set reset states back to what it was
Hd.PersistentMemory = resetval;
Hd.NumSamplesProcessed = Hd.Stage(1).NumSamplesProcessed;
