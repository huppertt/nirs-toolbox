function [isvalid, errmsg, errid] = thisvalidate(this)
%THISVALIDATE   Returns true if this object is valid.

%   Author(s): V. Pellissier
%   Copyright 2005 The MathWorks, Inc.

isvalid = 1;
errmsg = [];
errid = [];
if (this.NormalizedFrequency && this.Fpass~=1) || (~this.NormalizedFrequency && this.Fpass~=this.Fs/2)
    [isvalid, errmsg, errid] = checkincfreqs(this,{'Fpass','Fstop'});
end

% [EOF]
