function Hd2 = tolatticearma(Hd);
%TOLATTICEARMA  Convert to lattice ARMA.
%   Hd2 = TOLATTICEARMA(Hd) converts discrete-time filter Hd to lattice ARMA
%   filter Hd2.

%   Copyright 1988-2002 The MathWorks, Inc.
  
[b,a] = tf(Hd);
[k,v] = tf2latc(b,a);
Hd2 = dfilt.latticearma(k,v);