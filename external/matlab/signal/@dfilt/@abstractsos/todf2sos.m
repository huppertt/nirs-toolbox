function Hd2 = todf2sos(Hd);
%TODF2SOS  Convert to direct-form II sos.
%   Hd2 = TODF2SOS(Hd) converts discrete-time filter Hd to direct-form II
%   sos filter Hd2. 

%   Copyright 1988-2002 The MathWorks, Inc.
  

Hd2 = dfilt.df2sos(Hd.sosMatrix,Hd.ScaleValues);

