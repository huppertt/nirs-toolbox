function [isvalid, errmsg, errid] = thisvalidate(h)
%THISVALIDATE   

%   Author(s): R. Losada
%   Copyright 2003-2004 The MathWorks, Inc.

[isvalid, errmsg, errid] = checkincfreqs(h,{'Fstop','Fpass'});

% [EOF]
