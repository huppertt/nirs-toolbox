function designobj = getdesignobj(~, str)
%GETDESIGNOBJ Get the design object.

%   Copyright 2010 The MathWorks, Inc.

%#function fdfmethod.lpnormsbarbgrpdelay
designobj.iirlpnorm = 'fdfmethod.lpnormsbarbgrpdelay';

if nargin > 1
    designobj = designobj.(str);
end

% [EOF]
