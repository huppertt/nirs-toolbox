function [y,zf] = df2tfilter(q,b,a,x,zi)
% DF2TFILTER Filter for DFILT.DF2T class in single precision mode

%   Author(s): V.Pellissier
%   Copyright 1999-2004 The MathWorks, Inc.

x = quantizeinput(q,x);

[y,zf] = sdf2tfilter(b,a,x,zi);

