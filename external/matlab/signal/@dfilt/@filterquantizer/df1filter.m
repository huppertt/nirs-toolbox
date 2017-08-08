function [y,zfNum,zfDen,nBPtrf,dBPtrf] = df1filter(q,b,a,...
    x,ziNum,ziDen,nBPtr,dBPtr)
%DF1FILTER Filter for DFILT.DF1 class in double precision mode

%   Author(s): V.Pellissier
%   Copyright 1999-2004 The MathWorks, Inc.

x = quantizeinput(q,x);

% Call the DF1 filter implementation DLL
[y,zfNum,zfDen,nBPtrf,dBPtrf] = df1filter(b,a,x,ziNum,ziDen,nBPtr,dBPtr);


