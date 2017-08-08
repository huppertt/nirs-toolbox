%DOUBLE   Convert a FILTSTATES.DFIIR object to a double vector.
%   DOUBLE(Hs) returns the double precision value for the FILTSTATES.DFIIR
%   object.
%
%   EXAMPLE:
%   [b,a] = butter(4,.5); 
%   Hd = dfilt.df1(b,a);
%   hs = Hd.states; % Returns a FILTSTATES.DFIIR object
%   dblstates = double(hs); % Convert object to a double vector.
%
%   See also FILTSTATES/DFIIR, DFILT.

%   Author(s): P. Costa
%   Copyright 1988-2003 The MathWorks, Inc.

% [EOF]
