function [y,zf] = statespacefilter(q,Hd,x,zi)
% STATESPACEFILTER Filter for DFILT.STATESPACE class in double precision mode

%   Author(s): V.Pellissier
%   Copyright 1999-2004 The MathWorks, Inc.

x = quantizeinput(q,x);

if isempty(zi),
    y=Hd.D*x;
    zf = [];
    return
end

[y,zf] = statespacefilter(Hd.A,Hd.B,Hd.C,Hd.D,x,zi);

