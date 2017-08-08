function a = kratio(m,krat)
%KRATIO Utility function for use with ELLIP.
%   KRATIO(m,krat) is a function used to calculate the zeros of an
%   elliptic filter.  It is used with FMINSEARCH to find a parameter m 
%   satisfying ellipke(m)/ellipke(1-m) = krat.

%   Copyright 1988-2004 The MathWorks, Inc.

% to ensure we don't call ellipke(1) which is inf on non-ieee machines
% and that we only call with positive m.
m = min(1,max(m,0));
if abs(m) > eps && abs(m)+eps < 1
	k = ellipke([m,1-m]);
	r = k(1)./k(2) - krat;
elseif abs(m) <= eps	% m==0
	r = -krat;
else	% m==1 => r == inf, but can't for non-ieee machines
	r = 1e20;
end
a = abs(r);
