function cols = ncolors(n)

% Copyright 2011 The MathWorks, Inc.

d = 6/n;
hues = rem(10-d*(0:(n-1)),6);
cols = nncolor.hue2color(hues);
cols = bsxfun(@rdivide,cols,sum(cols,2));
