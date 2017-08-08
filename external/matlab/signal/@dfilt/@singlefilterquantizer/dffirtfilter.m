function [y,zf] = dffirtfilter(q,b,x,zi)

%   Author(s): R. Losada
%   Copyright 1988-2004 The MathWorks, Inc.


x = quantizeinput(q,x);
[y,zf] = sdffirtfilter(b,x,zi);

