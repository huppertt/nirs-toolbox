function Hd2 = todf1t(Hd);
%TODF1T  Convert to direct-form 1 transposed.
%   Hd2 = TODF1T(Hd) converts discrete-time filter Hd to direct-form 1
%   transposed filter Hd2.  

%   Copyright 1988-2002 The MathWorks, Inc.
  
[b,a] = tf(Hd);
Hd2 = dfilt.df1t(b,a);


