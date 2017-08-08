function Hd2 = todf2(Hd);
%TODF2  Convert to direct-form 2.
%   Hd2 = TODF2(Hd) converts discrete-time filter Hd to direct-form 2 filter
%   Hd2.  

%   Copyright 1988-2002 The MathWorks, Inc.
  
[b,a] = tf(Hd);
Hd2 = dfilt.df2(b,a);

