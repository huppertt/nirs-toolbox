function num = setrefnum(Hd, num)
%SETREFNUM Overloaded set on the refnum property.
  
%   Author: V. Pellissier
%   Copyright 1988-2003 The MathWorks, Inc.

validaterefcoeffs(Hd.filterquantizer, 'Numerator', num);
