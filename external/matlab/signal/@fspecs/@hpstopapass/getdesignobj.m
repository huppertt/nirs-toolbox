function designobj = getdesignobj(this, str)
%GETDESIGNOBJ   Get the design object.

%   Author(s): J. Schickler
%   Copyright 1988-2005 The MathWorks, Inc.

if isfdtbxinstalled
    %#function fdfmethod.eqriphpapass
    designobj.equiripple = 'fdfmethod.eqriphpapass';
else
    designobj = [];
end

if nargin > 1
    designobj = designobj.(str);
end

% [EOF]
