function k = local_max(x)
% k = local_max(x)
% finds location of local maxima
%
%   Copyright 1988-2005 The MathWorks, Inc.
s = size(x); x = [x(:)].'; N = length(x);
b1 = x(1:N-1)<=x(2:N); b2 = x(1:N-1)>x(2:N);
k = find(b1(1:N-2)&b2(2:N-1))+1;
if x(1)>x(2), k = [k, 1]; end
if x(N)>x(N-1), k = [k, N]; end
k = sort(k); if s(2) == 1, k = k'; end
