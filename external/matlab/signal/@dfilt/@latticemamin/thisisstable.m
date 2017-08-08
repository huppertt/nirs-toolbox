function f = thisisstable(Hd)
%THISISSTABLE  True if stable.
%   THISISSTABLE(Hd) returns 1 if discrete-time filter Hd is stable, and 0
%   otherwise. 
%
%   See also DFILT.   

%   Author(s): R. Losada
%   Copyright 1988-2003 The MathWorks, Inc.

% This should be private

f = true;  % nonrecursive filters are always stable.