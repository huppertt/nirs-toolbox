function designobj = getdesignobj(this, str)
%GETDESIGNOBJ   Get the designobj.

%   Copyright  The MathWorks, Inc.

if isfdtbxinstalled
    
    %#function fdfmethod.eqriphphbmin
    %#function fdfmethod.kaiserhphbmin
    %#function fdfmethod.elliphphalfbandmin
    %#function fdfmethod.butterhphalfbandmin
    %#function fdfmethod.iirhphalfbandeqripmin
    designobj.equiripple = 'fdfmethod.eqriphphbmin';
    designobj.kaiserwin  = 'fdfmethod.kaiserhphbmin';
    designobj.ellip      = 'fdfmethod.elliphphalfbandmin';
    designobj.butter     = 'fdfmethod.butterhphalfbandmin';
    designobj.iirlinphase = 'fdfmethod.iirhphalfbandeqripmin';
else
    designobj = [];
end

if nargin > 1
    designobj = designobj.(str);
end

% [EOF]
