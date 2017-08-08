function designobj = getdesignobj(this, str)
%GETDESIGNOBJ   Get the design object.

%   Copyright 2008 The MathWorks, Inc.

%#function fdfmethod.ellipparameqaudioshelfs

designobj.ellip = 'fdfmethod.ellipparameqaudioshelfs';

if nargin > 1
    designobj = designobj.(str);
end

% [EOF]
