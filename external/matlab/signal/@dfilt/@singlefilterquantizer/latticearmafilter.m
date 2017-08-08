function [y,zf] = latticearmafilter(q,k,kconj,ladder,x,zi)
% LATTICEARMAFILTER Filter for DFILT.LATTICEARMA class in single precision mode

%   Author(s): V.Pellissier
%   Copyright 1999-2004 The MathWorks, Inc.

x = quantizeinput(q,x);
[y,zf] = slatticearmafilter(k,kconj,ladder,x,zi);
