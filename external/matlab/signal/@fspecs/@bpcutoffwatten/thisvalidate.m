function [isvalid, errmsg, errid] = thisvalidate(h)
%THISVALIDATE   Check that this object is valid.

%   Copyright 2008 The MathWorks, Inc.

[isvalid, errmsg, errid] = checkincfreqs(h,{'Fcutoff1','Fcutoff2'});

% [EOF]
