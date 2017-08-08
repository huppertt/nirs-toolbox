function f = thisisstable(Hd)
%THISISSTABLE  True if filter is stable.
%   THISISSTABLE(Hd) returns 1 if discrete-time filter Hd is stable, and 0
%   otherwise. 
%
%   See also DFILT.   
  
%   Author: Thomas A. Bryan
%   Copyright 1988-2002 The MathWorks, Inc.

% This should be private

A = Hd.A;
if isempty(A)
  % If the state transition matrix is empty, then the filter is stable
  % because there are no states.
  f = true;
else
  f = max(abs(eig(A)))<1;
end
    
           

