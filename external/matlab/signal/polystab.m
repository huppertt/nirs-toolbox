function b = polystab(a)
%POLYSTAB Polynomial stabilization.
%   POLYSTAB(A), where A is a vector of polynomial coefficients,
%   stabilizes the polynomial with respect to the unit circle;
%   roots whose magnitudes are greater than one are reflected
%   inside the unit circle.
%
%   % Example:
%   %   Convert a linear-phase filter into a minimum-phase filter with the 
%   %   same magnitude response.
%
%   h = fir1(25,0.4);               % Window-based FIR filter design
%   flag_linphase = islinphase(h)   % Determines if filter is linear phase
%   hmin = polystab(h) * norm(h)/norm(polystab(h)); 
%   flag_minphase = isminphase(hmin)% Determines if filter is minimum phase

%   Author(s): J.N. Little,7-25-89, handles roots at zero
%   Copyright 1988-2004 The MathWorks, Inc.

if isempty(a), b = a; return, end
if length(a) == 1, b = a; return, end
v = roots(a); i = find(v~=0);
vs = 0.5*(sign(abs(v(i))-1)+1);
v(i) = (1-vs).*v(i) + vs./conj(v(i));
ind = find(a~=0);
b = a(ind(1))*poly(v);

% Return only real coefficients if input was real:
if ~any(imag(a))
	b = real(b);
end

