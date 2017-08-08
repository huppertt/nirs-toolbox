function den = getdenominator(Hd, den)
%GETDENOMINATOR Overloaded get on the Denominator property.
  
%   Copyright 1988-2003 The MathWorks, Inc.

den = getdenominator(Hd.filterquantizer,Hd.privden);
