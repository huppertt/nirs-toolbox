function [y,zf] = latticeallpassfilter(q,k,kconj,x,zi)
% LATTICEALLPASSFILTER Filter for DFILT.LATTICEALLPASS class in double precision mode

%   Author(s): V.Pellissier
%   Copyright 1999-2004 The MathWorks, Inc.

x = quantizeinput(q,x);
[y,zf] = latticeallpassfilter(k,kconj,x,zi);

