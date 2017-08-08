function [y,zf] = secfilter(Hd,x,zi)
%SECFILTER Filter this section.
%   [Y,Zf] = SECFILTER(Hd,X,ZI) filters this section.
%
%   See also DFILT.    

%   Author(s): R. Losada
%   Copyright 2005 The MathWorks, Inc.

  
q = Hd.filterquantizer;
a = Hd.privallpasscoeffs;
c = wdfcoefficients(Hd);
[y,zf] = wdfallpassfilter(q,c.Section1,x,zi);



% [EOF]
