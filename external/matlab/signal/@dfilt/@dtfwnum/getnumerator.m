function num = getnumerator(Hd, num)
%GETNUMERATOR Overloaded get on the Numerator property.
  
%   Copyright 1988-2003 The MathWorks, Inc.

num = getnumerator(Hd.filterquantizer, Hd.privnum);
