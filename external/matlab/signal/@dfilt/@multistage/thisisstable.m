function f = isstable(Hd)
%ISSTABLE  True if filter is stable.
%   ISSTABLE(Hd) returns 1 if discrete-time filter Hd is stable, and 0
%   otherwise. 
%
%   See also DFILT.   

%   Author: Thomas A. Bryan, J. Schickler
%   Copyright 1988-2004 The MathWorks, Inc.

% This should be private

f = all(isstable(Hd.Stage));
