function designobj = getdesignobj(~, str)
%GETDESIGNOBJ   Get the design object.

%   Copyright 2010 The MathWorks, Inc.

%#function fdfmethod.lpnormmultibandarbgrpdelay
designobj.iirlpnorm = 'fdfmethod.lpnormmultibandarbgrpdelay';

if nargin > 1
    designobj = designobj.(str);
end

% [EOF]
