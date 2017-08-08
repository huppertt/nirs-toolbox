function validatedata(this, data)
%VALIDATEDATA   Validate the data for this object.

%   Author(s): P. Pacheco
%   Copyright 1988-2006 The MathWorks, Inc.

% Call "private" function (using a trick) which is also used by abstractps.
dspdata.validatedata(this,data);

% [EOF]
