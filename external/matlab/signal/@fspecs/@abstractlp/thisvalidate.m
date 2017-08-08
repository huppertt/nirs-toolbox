function [isvalid, errmsg, errid] = thisvalidate(h)
%THISVALIDATE   Validate this object.

%   Author(s): R. Losada
%   Copyright 2003-2005 The MathWorks, Inc.

[isvalid, errmsg, errid] = checkincfreqs(h, {'Fpass','Fstop'});

% [EOF]
