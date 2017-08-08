function f = thisisstable(Hd)
%THISISSTABLE  True if filter is stable.
%   THISISSTABLE(Hd) returns 1 if discrete-time filter Hd is stable, and 0
%   otherwise. 

%   Author(s): R. Losada, T. Bryan
%   Copyright 1988-2002 The MathWorks, Inc.

% This should be private

K = Hd.Lattice;
if isempty(K)
  % If the lattice is empty, then the filter is stable because there is no
  % feedback. 
  f = true;
else
  % Nonempty Lattice AR and ARMA filters are stable iff the maximum magnitude
  % lattice coefficient is less than 1.
  f = max(abs(K))<1;
end          

