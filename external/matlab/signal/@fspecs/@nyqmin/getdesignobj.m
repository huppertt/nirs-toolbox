function designobj = getdesignobj(this, str)
%GETDESIGNOBJ   Get the designobj.

%   Author(s): J. Schickler
%   Copyright 1988-2005 The MathWorks, Inc.

if isfdtbxinstalled
    
    %#function fdfmethod.kaisernyqmin
    %#function fdfmethod.eqripnyqmin
    %#function fdfmethod.multistagenyq
    designobj.kaiserwin  = 'fdfmethod.kaisernyqmin';
    designobj.equiripple = 'fdfmethod.eqripnyqmin';
    designobj.multistage = 'fdfmethod.multistagenyq';
else
    designobj = [];
end

if nargin > 1
    designobj = designobj.(str);
end

% [EOF]
