function v = vratio(u,ineps,mp)
%VRATIO Utility function for use with ELLIP.
%   VRATIO(u,ineps,mp) is a function used to calculate the poles of an
%   elliptic filter.  It finds a u so sn(u)/cn(u) = 1/epsilon ( = ineps), 
%   with parameter mp.

%   Copyright 1988-2002 The MathWorks, Inc.

%   global information - 1/epsilon, the value s/c should attain
%   with parameter mp.

[s,c] = ellipj(u,mp);
v = abs(ineps - s/c);

