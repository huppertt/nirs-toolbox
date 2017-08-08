function C = setc(Hd, C)
%SETC Overloaded set function on the C property.
  
%   Author: V. Pellissier
%   Copyright 1988-2003 The MathWorks, Inc.
  
set(Hd,'refC',C);
quantizecoeffs(Hd);

% Hold an empty to not duplicate storage
C = [];

