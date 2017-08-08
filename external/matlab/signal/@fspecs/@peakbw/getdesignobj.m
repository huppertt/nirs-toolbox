function designobj = getdesignobj(this, str)
%GETDESIGNOBJ   Get the design object.

%   Author(s): J. Schickler
%   Copyright 2006 The MathWorks, Inc.

%#function fdfmethod.butterpeakbw

designobj.butter = 'fdfmethod.butterpeakbw';

if nargin > 1
    designobj = designobj.(str);
end

% [EOF]
