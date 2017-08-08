function specs = whichspecs(h)
%WHICHSPECS

%   Author(s): J. Schickler
%   Copyright 1988-2002 The MathWorks, Inc.

specs = ft_whichspecs(h);

newspecs.name     = 'invSincFreqFactor';
newspecs.datatype = 'udouble';
newspecs.defval   = 1;
newspecs.callback = [];
newspecs.descript = 'magspec';

specs = [newspecs specs];

% [EOF]
