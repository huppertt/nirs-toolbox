function designobj = getdesignobj(this, str)
%GETDESIGNOBJ   Get the design object.

%   Author(s): V. Pellissier
%   Copyright 2005 The MathWorks, Inc.

%#function fdfmethod.lpnormbp2
designobj.iirlpnorm = 'fdfmethod.lpnormbp2';

if nargin > 1
    designobj = designobj.(str);
end

% [EOF]
