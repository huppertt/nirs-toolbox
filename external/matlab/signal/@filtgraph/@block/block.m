function B = block(index,blktype,lbl)
%BLOCK Constructor for this class.

%   Author(s): Roshan R Rammohan
%   Copyright 1988-2005 The MathWorks, Inc.

error(nargchk(0,3,nargin,'struct'));

B = filtgraph.block;

if nargin > 0
    B.nodeIndex = index;
end

if nargin > 1
    B.blocktype = blktype;
    B.label = blktype;

    if ~(strcmpi(blktype,'output')||strcmpi(blktype,'goto'))
        B.setnumoutports(1);
    end

    if ~(strcmpi(blktype,'input')||strcmpi(blktype,'from'))
        if ~(strcmpi(blktype,'sum') || strcmpi(blktype,'mult'))
            B.setnuminports(1);
        else
            B.setnuminports(2);
        end
    end

end

if nargin > 2
    B.label = lbl;
end

