function [isvalid, errmsg, errid] = thisvalidate(h)
%THISVALIDATE Check that this object is valid.

%   Copyright 2011 The MathWorks, Inc.

[isvalid, errmsg, errid] = checkincfreqs(h,{'Fstop','Fpass'});

% [EOF]
