function newval = todb(h,val,notused)
%TODB Convert squared value to dB.

%   Author(s): R. Losada
%   Copyright 1988-2002 The MathWorks, Inc.

newval = 10*log10(1/val);
