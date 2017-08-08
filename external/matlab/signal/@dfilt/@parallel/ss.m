function [A,B,C,D] = ss(Hd)
  
%SS  Discrete-time filter to state-space conversion.
%   [A,B,C,D] = SS(Hd) converts discrete-time filter Hd to state-space
%   representation given by 
%     x(k+1) = A*x(k) + B*u(k)
%     y(k)   = C*x(k) + D*u(k)
%   where x is the state vector, u is the input vector, and y is the output
%   vector. 
%
%   See also DFILT.

%   Author(s): R. Losada, T. Bryan
%   Copyright 1988-2006 The MathWorks, Inc.

% Check if there are multirate filters in the parallel
for k = 1:nstages(Hd),
    if any(strcmpi(class(Hd.Stage(k)),{'mfilt.abstractmultirate','mfilt.cascade'})),
        error(message('signal:dfilt:parallel:ss:invalidFilters'));
    end
end

% Initialize state-space model
[A,B,C,D] = ss(Hd.Stage(1));

% Form the state-space model of the parallel recursively
for k = 2:length(Hd.Stage)
    
  % Generate state-space model per section
  [a2,b2,c2,d2] = ss(Hd.Stage(k));
   
  % Combine section with overall state-space model
  A = [A,                           zeros(size(A,1),size(a2,2));
       zeros(size(a2,1),size(A,2)),        a2];
  % Force uniformity with empty state matrix.
  if isempty(A)
    A = [];
    B = zeros(0,1);
    C = zeros(1,0);
  else
    B = [B;b2];
    C = [C, c2];
  end
  D = D + d2;
end

