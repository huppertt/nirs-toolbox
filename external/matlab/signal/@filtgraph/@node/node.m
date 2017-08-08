function N = node(blki,plist)
%NODE Constructor for this class.
% If blki is not of type filtgraph.block then a block of type blki can be
% created

%   Author(s): Roshan R Rammohan
%   Copyright 1988-2004 The MathWorks, Inc.

error(nargchk(0,2,nargin,'struct'));

N = filtgraph.node;

if nargin > 0
    if ~isa(blki,'filtgraph.block')
        blki = filtgraph.block(N.index,blki);
    end
    blki = copy(blki);
    N.block = blki;
end

if nargin > 1
    N.paramList = plist;
end
