function B = setb(Hd, B)
%SETB Overloaded set function on the B property.
  
%   Author: V. Pellissier
%   Copyright 1988-2003 The MathWorks, Inc.
  
set(Hd,'refB',B);
quantizecoeffs(Hd);

% Hold an empty to not duplicate storage
B = [];
