function designobj = getdesignobj(this, str)
%GETDESIGNOBJ   Get the design object.

%   Author(s): V. Pellissier
%   Copyright 2006 The MathWorks, Inc.

%#function fdfmethod.butteroctave
designobj.butter = 'fdfmethod.butteroctave';

if nargin > 1
    designobj = designobj.(str);
end

% [EOF]
