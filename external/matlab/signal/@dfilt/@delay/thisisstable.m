function f = thisisstable(Hd)
%THISISSTABLE  True if filter is stable.
%   THISISSTABLE(Hd) returns 1 if discrete-time filter Hd is stable, and 0
%   otherwise. 

%   Author(s): T. Bryan
%   Copyright 1988-2005 The MathWorks, Inc.

% This should be private

% Scalar filters are always stable.
f = true;
            

