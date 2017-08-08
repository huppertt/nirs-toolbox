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
%   Copyright 1988-2004 The MathWorks, Inc.


if isempty(Hd.Stage),
  error(message('signal:dfilt:cascade:ss:DFILTErr'));
end

% Form the state-space model of the cascade recursively
[A,B,C,D] = ss(Hd.Stage(1));

% Turn warnings off in case empty matrices are produced, which are okay, but
% will warn if an empty is multiplied by a non-empty.
wrn = warning('off');
for k = 2:length(Hd.Stage)
  % Generate state-space model per section
  [a2,b2,c2,d2] = ss(Hd.Stage(k));
   
  % Combine section with overall state-space model
  A = [A, zeros(size(A,1),size(a2,2));b2*C, a2];
  if isempty(A)
    A = [];
    B = zeros(0,1);
    C = zeros(1,0);
  else
    B = [B;b2*D];
    C = [d2*C, c2];
  end
  D = d2*D;
end
% Restore the warning state
warning(wrn);
