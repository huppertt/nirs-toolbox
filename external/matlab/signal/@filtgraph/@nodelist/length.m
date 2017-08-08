function L = length(NodeList)
%LENGTH of NodeList

%   Author(s): Roshan R Rammohan
%   Copyright 1988-2004 The MathWorks, Inc.

error(nargchk(1,1,nargin,'struct'));

NL = NodeList;
L = length(NL.nodes);
