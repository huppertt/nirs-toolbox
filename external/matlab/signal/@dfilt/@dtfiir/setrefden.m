function den = setrefden(Hd, den)
%SETREFNUM Overloaded set on the refden property.
  
%   Author: P. Costa
%   Copyright 1988-2004 The MathWorks, Inc.

validaterefcoeffs(Hd.filterquantizer, 'Denominator', den);
