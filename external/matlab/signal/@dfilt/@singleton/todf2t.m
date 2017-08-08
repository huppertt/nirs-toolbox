function Hd2 = todf2t(Hd);
%TODF2T  Convert to direct-form 2 transposed.
%   Hd2 = TODF2T(Hd) converts discrete-time filter Hd to direct-form 2
%   transposed filter Hd2.  

%   Copyright 1988-2002 The MathWorks, Inc.
  
[b,a] = tf(Hd);
Hd2 = dfilt.df2t(b,a);