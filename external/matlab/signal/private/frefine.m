function rs = frefine(a,v,rs)
% f = frefine(a,v,rs);
% refine local minima and maxima of H using Newton's method
%
% H  : H = cos(w*v)*a;
% rs : initial values for the extrema of H
%
%   Copyright 1988-2002 The MathWorks, Inc.
a = a(:);
v = v(:)';
w = rs(:);
m = length(a)-1;
for k = 1:5
   H1 = -sin(w*v) * (v'.*a);
   H2 = -cos(w*v) * ((v.^2)'.*a);
   w = w - H1./H2;
end
rs(:) = w;
