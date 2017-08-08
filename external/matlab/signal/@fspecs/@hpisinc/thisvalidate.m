function [isvalid, errmsg, errid] = thisvalidate(h)
%THISVALIDATE   

%   Copyright 2011 The MathWorks, Inc.

[isvalid, errmsg, errid] = checkincfreqs(h,{'Fstop','Fpass'});

% [EOF]
