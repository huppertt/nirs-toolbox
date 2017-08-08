function [y,zf] = secfilter(this,x,zi)
%SECFILTER Filter this section.
%   [Y,Zf] = SECFILTER(this,X,ZI) filters this section.  This function is only
%   intended to be called from DFILT/FILTER.  

%   Author(s): V. Pellissier, M.Chugh
%   Copyright 2005 The MathWorks, Inc.

q = this.filterquantizer;
b = this.Latency;
[y,zf] = delayfilter(q,b,x,zi);

% [EOF]
