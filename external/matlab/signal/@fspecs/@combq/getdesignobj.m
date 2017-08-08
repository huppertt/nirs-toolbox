function designobj = getdesignobj(this, str)
%GETDESIGNOBJ   Get the design object.

%   Copyright 2008 The MathWorks, Inc.

%#function fdfmethod.buttercombq

designobj.butter = 'fdfmethod.buttercombq';

if nargin > 1
    designobj = designobj.(str);
end

% [EOF]
