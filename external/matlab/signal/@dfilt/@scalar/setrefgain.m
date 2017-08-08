function g = setrefgain(Hd, g)
%SETREFGAIN Overloaded set on the refgain property.
  
%   Author: V. Pellissier
%   Copyright 1988-2004 The MathWorks, Inc.

validaterefcoeffs(Hd.filterquantizer, 'Gain', g);
