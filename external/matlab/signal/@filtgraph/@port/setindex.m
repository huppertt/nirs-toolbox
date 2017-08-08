function P = setindex(Pi,index)

%   Author(s): Roshan R Rammohan
%   Copyright 1988-2004 The MathWorks, Inc.

error(nargchk(2,2,nargin,'struct'));
P = Pi;

P.nodeIndex = index;
