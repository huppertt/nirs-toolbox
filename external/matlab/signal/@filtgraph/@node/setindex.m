function N = setindex(Ni,index)
%SETINDEX

%   Author(s): Roshan R Rammohan
%   Copyright 1988-2004 The MathWorks, Inc.

error(nargchk(2,2,nargin,'struct'));

N = Ni;
N.index = index;

if ~(N.block.nodeIndex == index)
    N.block.setindex(index);
end
