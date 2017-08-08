function y = filter(Hd,x,dim)
%FILTER Discrete-time filter.
%   Y= FILTER(Hd,X,DIM) filters data X over dimension DIM with and returns
%   output Y. 
%
%   See also DFILT.   
  
%   Author: Thomas A. Bryan
%   Copyright 1988-2010 The MathWorks, Inc.
  
if nargin<3, dim=[]; end

if isempty(x), 
  y = x;
  return; 
end

% Because this filter method used to allow for specifying the states as the
% third input, we check that the DIM input is not a vector/matrix.
if any(size(dim)>1),
    error(message('signal:dfilt:cascade:filter:DimMustBeInt'));
end

% Get current value of reset states
resetval = Hd.PersistentMemory;

y = x;

flagdiff = false;
for k=1:length(Hd.Stage)
    if Hd.Stage(k).PersistentMemory~=resetval,
        Hd.Stage(k).PersistentMemory = resetval;
        flagdiff = true;
    end
    y = filter(Hd.Stage(k),y,dim);
end

if flagdiff,
     warning(message('signal:dfilt:cascade:filter:flagInconsistency', mat2str( resetval )));
end

% Set reset states back to what it was
Hd.PersistentMemory = resetval;

Hd.NumSamplesProcessed = Hd.Stage(1).NumSamplesProcessed;
