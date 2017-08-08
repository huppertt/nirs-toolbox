function [y,zf] = latticemaminfilter(q,k,kconj,x,zi)
% LATTICEMAMINFILTER Filter for DFILT.LATTICEMAMIN class in double precision mode

%   Author(s): V.Pellissier
%   Copyright 1999-2004 The MathWorks, Inc.

x = quantizeinput(q,x);
[y,zf] = latticemaminphasefilter(k,kconj,x,zi);

