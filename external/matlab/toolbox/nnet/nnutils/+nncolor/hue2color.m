function col = hue2color(hue)

% Copyright 2011 The MathWorks, Inc.

r = max(0,min(1,abs(hue-3)-1));
hue = rem(hue+4,6);
g = max(0,min(1,abs(hue-3)-1));
hue = rem(hue+4,6);
b = max(0,min(1,abs(hue-3)-1));
col = [r' g' b'];
