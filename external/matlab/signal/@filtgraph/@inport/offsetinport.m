function iprt = offsetinport(iprt,offset)
%OFFSETINPORT Offset the from.node of an inport 

%   Author(s): S Dhoorjaty
%   Copyright 1988-2004 The MathWorks, Inc.


iprt.from = offsetnodeport(iprt.from,offset);
