function [y,zf] = latticearfilter(q,k,kconj,ladder,x,zi)
% LATTICEARMAFILTER Filter for DFILT.LATTICEARMA class in double precision mode

%   Author(s): V.Pellissier
%   Copyright 1999-2004 The MathWorks, Inc.

x = quantizeinput(q,x);
[y,zf] = latticearmafilter(k,kconj,ladder,x,zi);

