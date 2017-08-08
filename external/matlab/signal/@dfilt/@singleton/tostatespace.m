function Hd2 = tostatespace(Hd);
%TOSTATESPACE  Convert to statespace.
%   Hd2 = TOSTATESPACE(Hd) converts discrete-time filter Hd to statespace
%   filter Hd2.

%   Copyright 1988-2002 The MathWorks, Inc.
  
[A,B,C,D] = ss(Hd);
Hd2 = dfilt.statespace(A,B,C,D);