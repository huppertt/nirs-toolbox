function designobj = getdesignobj(this, str)
%GETDESIGNOBJ   Get the design object.

%   Author(s): J. Schickler
%   Copyright 1988-2006 The MathWorks, Inc.

%#function fmethod.cheby2bs
designobj.cheby2 = 'fmethod.cheby2bs';

if isfdtbxinstalled
    %#function fdfmethod.ellipbscutoffwas
    designobj.ellip  = 'fdfmethod.ellipbscutoffwas';
end

if nargin > 1
    designobj = designobj.(str);
end

% [EOF]
