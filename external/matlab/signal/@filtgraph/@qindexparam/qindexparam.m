function indparm = qindexparam(index,paramlist)
%ASSOC Constructor for this class.

%   Author(s): Roshan R Rammohan
%   Copyright 1988-2004 The MathWorks, Inc.

error(nargchk(2,2,nargin,'struct'));

indparm = filtgraph.qindexparam;

indparm.index = index;

indparm.params = paramlist;
