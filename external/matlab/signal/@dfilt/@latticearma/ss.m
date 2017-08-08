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
%   Copyright 1988-2002 The MathWorks, Inc.


k = Hd.Lattice;
v = Hd.Ladder;

% Handle the scalar case
if length(k)==1 & length(v) == 1 & k == 0,
    A = [];
    B = zeros(0,1);
    C = zeros(1,0);
    D = v;
    return
end

k = k(:);
v = v(:);

% length(v) must be length(k)+1.  Fill with zeros at the end if necessary.
N = length(k);
M = length(v);
if M ~= N+1
  % They are not the right length
  if M<=N
    % v is too short
    v = [v;zeros(N-M+1,1)];
    M = length(v);
  else
    % k is too short
    k = [k;zeros(M-N,1)];
    N = length(k);
  end
end

[A,B] = statespaceab(Hd);

% Force uniformity with empty state matrix.
if isempty(A)
  A = [];
  B = zeros(0,1);
  C = zeros(1,0);
else
  C = sum(v(1:N,ones(N,1)).*A,1);
  C(end) = C(end) + (v(end)*(1-abs(k(end)).^2));
end

D = sum([v(1);v(2:end).*conj(k)]);
