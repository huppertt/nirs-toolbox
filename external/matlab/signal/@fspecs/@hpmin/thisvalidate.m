function [isvalid, errmsg, errid] = thisvalidate(h)
%THISVALIDATE   Check that this object is valid.

%   Author(s): R. Losada
%   Copyright 1988-2004 The MathWorks, Inc.

[isvalid, errmsg, errid] = checkincfreqs(h,{'Fstop','Fpass'});

% [EOF]
