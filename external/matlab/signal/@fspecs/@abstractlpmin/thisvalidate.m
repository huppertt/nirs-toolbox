function [isvalid, errmsg, errid] = thisvalidate(h)
%THISVALIDATE   Checks if this object is valid.

%   Author(s): R. Losada
%   Copyright 1988-2005 The MathWorks, Inc.

[isvalid, errmsg, errid] = checkincfreqs(h,{'Fpass','Fstop'});

% [EOF]
