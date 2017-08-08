function designobj = getdesignobj(this, str)
%GETDESIGNOBJ   Get the designobj.

%   Copyright 2007 The MathWorks, Inc.

if isfdtbxinstalled
    
    %#function fdfmethod.eqriphphbordntw
    %#function fdfmethod.kaiserhphbordntw
    %#function fdfmethod.firlshphbordntw
    %#function fdfmethod.elliphphalfbandfpass
    %#function fdfmethod.iirhphalfbandeqripfpass
    designobj.equiripple  = 'fdfmethod.eqriphphbordntw';
    designobj.kaiserwin   = 'fdfmethod.kaiserhphbordntw';
    designobj.firls       = 'fdfmethod.firlshphbordntw';
    designobj.ellip       = 'fdfmethod.elliphphalfbandfpass';
    designobj.iirlinphase = 'fdfmethod.iirhphalfbandeqripfpass';
else
    designobj = [];
end

if nargin > 1
    designobj = designobj.(str);
end

% [EOF]
