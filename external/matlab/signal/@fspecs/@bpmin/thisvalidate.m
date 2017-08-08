function [isvalid, errmsg, errid] = thisvalidate(h)
%THISVALIDATE   Checks if this object is valid.

%   Author(s): R. Losada
%   Copyright 1988-2004 The MathWorks, Inc.

[isvalid, errmsg, errid] = checkincfreqs(h,{'Fstop1','Fpass1','Fpass2','Fstop2'});

% [EOF]
