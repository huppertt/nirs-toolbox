function designobj = getdesignobj(this, str)
%GETDESIGNOBJ   Get the design object.

%   Author(s): V. Pellissier
%   Copyright 2005-2006 The MathWorks, Inc.

%#function fdfmethod.ellipbscutoffwapas
designobj.ellip = 'fdfmethod.ellipbscutoffwapas';

if nargin > 1
    designobj = designobj.(str);
end

% [EOF]
