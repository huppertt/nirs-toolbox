function [y,zf,tapIndex] = dfantisymmetricfirfilter(q,b,x,zi,tapIndex)
% DFANTISYMMETRICFIRFILTER Filter for DFILT.DFASYMFIR class in single precision mode

%   Author(s): V.Pellissier
%   Copyright 1999-2004 The MathWorks, Inc.

x = quantizeinput(q,x);
[y,zf,tapIndex] = sdfantisymmetricfirfilter(b,x,zi,tapIndex);

