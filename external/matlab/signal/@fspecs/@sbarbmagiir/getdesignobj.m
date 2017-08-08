function designobj = getdesignobj(this, str)
%GETDESIGNOBJ   Get the design object.

%   Author(s): V. Pellissier
%   Copyright 2005 The MathWorks, Inc.

%#function fdfmethod.lpnormsbarbmag2
designobj.iirlpnorm = 'fdfmethod.lpnormsbarbmag2';

if nargin > 1
    designobj = designobj.(str);
end

% [EOF]
