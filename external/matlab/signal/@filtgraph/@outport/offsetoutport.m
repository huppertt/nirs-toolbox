function oprt = offsetoutport(oprt,n,offset)
%OFFSETOUTPORT Offset the to(n).node of outport oprt

%   Author(s): S Dhoorjaty
%   Copyright 1988-2004 The MathWorks, Inc.


oprt.to(n) = offsetnodeport(oprt.to(n),offset);
