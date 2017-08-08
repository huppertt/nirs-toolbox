function Hd2 = todf1(Hd);
%TODF1  Convert to direct-form 1.
%   Hd2 = TODF1(Hd) converts discrete-time filter Hd to direct-form 1 filter
%   Hd2. 

%   Copyright 1988-2002 The MathWorks, Inc.
  
[b,a] = tf(Hd);
Hd2 = dfilt.df1(b,a);

