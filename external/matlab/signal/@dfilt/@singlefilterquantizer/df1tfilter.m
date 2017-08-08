function [y,zfNum,zfDen] = df1tfilter(q,b,a,x,ziNum,ziDen)
% DF1TFILTER Filter for DFILT.DF1T class in single precision mode

%   Author(s): V.Pellissier
%   Copyright 1999-2004 The MathWorks, Inc.

x = quantizeinput(q,x);

[y,zfNum,zfDen] = sdf1tfilter(b,a,x,ziNum,ziDen);
