function [isvalid, errmsg, errid] = thisvalidate(h)
%THISVALIDATE   Check that this object is valid.

%   Author(s): R. Losada
%   Copyright 2003-2005 The MathWorks, Inc.

[isvalid, errmsg, errid] = checkincfreqs(h,{'Fpass','Fstop'});

% [EOF]
