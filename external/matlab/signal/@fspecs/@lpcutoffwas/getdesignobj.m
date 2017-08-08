function designobj = getdesignobj(this, str)
%GETDESIGNOBJ   Get the design object.

%   Author(s): J. Schickler
%   Copyright 1988-2006 The MathWorks, Inc.

%#function fmethod.cheby2lp
designobj.cheby2 = 'fmethod.cheby2lp';

if isfdtbxinstalled
    %#function fdfmethod.elliplpcutoffwas
    designobj.ellip  = 'fdfmethod.elliplpcutoffwas';
end

if nargin > 1
    designobj = designobj.(str);
end

% [EOF]
