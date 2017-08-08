function [y,zf,tapIndex] = dffirfilter(q,b,x,zi,tapIndex)

%   Author(s): R. Losada
%   Copyright 1988-2004 The MathWorks, Inc.

x = quantizeinput(q,x);
[y,zf,tapIndex] = dffirfilter(b,x,zi,tapIndex);
