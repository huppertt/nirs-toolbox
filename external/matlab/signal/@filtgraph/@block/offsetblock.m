function blk = offsetblock(blk,offset)
%OFFSETBLOCK Offset the nodeIndex of block blk and its contents

%   Author(s): S Dhoorjaty
%   Copyright 1988-2005 The MathWorks, Inc.

blk.nodeIndex = blk.nodeIndex + offset;

if ~isempty(blk.outport)
    
    for outpm = 1:length(blk.outport)
        blk.outport(outpm) = offsetport(blk.outport(outpm), offset);
    end

end

for m = 1:length(blk.inport)
    blk.inport(m) = offsetport(blk.inport(m),offset);
end
