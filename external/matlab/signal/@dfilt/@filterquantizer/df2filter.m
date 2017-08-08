function [y,zf,tapidxf] = df2filter(q,b,a,x,zi,tapidxi)
% DF2FILTER Filter for DFILT.DF2 class in double precision mode

%   Author(s): V.Pellissier
%   Copyright 1999-2004 The MathWorks, Inc.

x = quantizeinput(q,x);

% Call the DF2 filter implementation DLL
[y,zf,tapidxf] = df2filter(b,a,x,zi,tapidxi);
