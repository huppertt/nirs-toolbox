function [isvalid, errmsg, errid] = thisvalidate(h)
%THISVALIDATE   Checks if this object is valid.

%   Author(s): R. Losada
%   Copyright 1988-2004 The MathWorks, Inc.

[isvalid, errmsg, errid] = checkincfreqs(h,{'Fpass1','Fstop1','Fstop2','Fpass2'});

% [EOF]
