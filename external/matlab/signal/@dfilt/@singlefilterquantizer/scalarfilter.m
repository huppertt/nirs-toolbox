function [y,zf] = scalarfilter(q,b,x,zi)
% SCALARFILTER Filter for DFILT.SCALAR class in single precision mode

%   Author(s): V.Pellissier
%   Copyright 1999-2004 The MathWorks, Inc.

x = quantizeinput(q,x);
y = b * x;
zf = single(zi);

