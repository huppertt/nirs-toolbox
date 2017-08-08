function prt = offsetport(prt,offset)
%OFFSETPORT Offsets the nodeIndex of port prt

%   Author(s): S Dhoorjaty
%   Copyright 1988-2004 The MathWorks, Inc.

prt.nodeIndex = prt.nodeIndex + offset;
