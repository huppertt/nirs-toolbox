function designobj = getdesignobj(this, str)
%GETDESIGNOBJ   Get the design object.

%   Copyright 2008 The MathWorks, Inc.

%#function fdfmethod.buttercombbw

designobj.butter = 'fdfmethod.buttercombbw';

if nargin > 1
    designobj = designobj.(str);
end

% [EOF]
