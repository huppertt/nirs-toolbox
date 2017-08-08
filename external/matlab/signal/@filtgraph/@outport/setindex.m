function op = setindex(OPi,index)
%SETINDEX

%   Author(s): Roshan R Rammohan
%   Copyright 1988-2004 The MathWorks, Inc.

error(nargchk(2,2,nargin,'struct'));

op = OPi;

op.nodeIndex = index;




