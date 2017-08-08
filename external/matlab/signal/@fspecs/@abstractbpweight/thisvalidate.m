function [isvalid, errmsg, errid] = thisvalidate(h)
%THISVALIDATE   

%   Copyright 2003-2011 The MathWorks, Inc.

[isvalid, errmsg, errid] = checkincfreqs(h,{'Fstop1','Fpass1','Fpass2','Fstop2'});

% [EOF]
