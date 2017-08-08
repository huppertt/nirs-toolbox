function s = this_setstage(Hd,s) 
%THIS_SETSTAGE PreSet function for the stage property.

%   Copyright 2009 The MathWorks, Inc.

% Create a listener for the stage
l  = handle.listener(s, 'clearmetadata', @clearmetadata); 
set(l,  'callbacktarget', Hd);
set(Hd, 'clearmetadatalistener', l);

clearmetadata(Hd);

