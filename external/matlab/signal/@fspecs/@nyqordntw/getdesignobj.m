function designobj = getdesignobj(this, str)
%GETDESIGNOBJ   Get the designobj.

%   Author(s): J. Schickler
%   Copyright 1988-2005 The MathWorks, Inc.

if isfdtbxinstalled        
    %#function fdfmethod.eqripnyqordntw
    %#function fdfmethod.kaiserhbordntw
    designobj.equiripple = 'fdfmethod.eqripnyqordntw';
    designobj.kaiserwin  = 'fdfmethod.kaiserhbordntw';
else
    designobj = [];
end

if nargin > 1
    designobj = designobj.(str);
end

% [EOF]
