function [isvalid, errmsg, errid] = thisvalidate(h)
%THISVALIDATE   

%   Copyright 2003-2011 The MathWorks, Inc.

[isvalid, errmsg, errid] = checkincfreqs(h,{'Fpass1','Fstop1','Fstop2','Fpass2'});

% [EOF]
