function Hd2 = todf1sos(Hd);
%TODF1SOS  Convert to direct-form 1 sos.
%   Hd2 = TODF1SOS(Hd) converts discrete-time filter Hd to direct-form 1
%   sos filter Hd2. 

%   Copyright 1988-2002 The MathWorks, Inc.
  

Hd2 = dfilt.df1sos(Hd.sosMatrix,Hd.ScaleValues);

