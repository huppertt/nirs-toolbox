function isstableflag = thisisstable(Hd)
%THISISSTABLE  True if filter is stable.
%   THISISSTABLE(Hd) returns 1 if discrete-time filter Hd is stable, and 0
%   otherwise. 

%   Copyright 1988-2012 The MathWorks, Inc.

% This should be private

isstableflag = isstable(Hd.Numerator,Hd.Denominator);

    
            

