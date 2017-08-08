function [y,zf] = latticemamaxfilter(q,k,kconj,x,zi)
% LATTICEMAMAXFILTER Filter for DFILT.LATTICEMAMAX class in double precision mode

%   Author(s): V.Pellissier
%   Copyright 1999-2004 The MathWorks, Inc.

x = quantizeinput(q,x);
[y,zf] = latticemamaxphasefilter(k,kconj,x,zi);

