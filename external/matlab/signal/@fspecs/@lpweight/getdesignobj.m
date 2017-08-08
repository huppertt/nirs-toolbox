function designobj = getdesignobj(~, str, sigonlyflag)
%GETDESIGNOBJ   Get the designobj.

%   Copyright 1988-2013 The MathWorks, Inc.

if nargin < 2
    str = [];
    sigonlyflag = false;
elseif nargin < 3
    sigonlyflag = false;
end

if isfdtbxinstalled && ~sigonlyflag
    %#function fdfmethod.eqriplp
    %#function fdfmethod.lpnormlp1
    designobj.equiripple = 'fdfmethod.eqriplp';
    designobj.iirlpnorm  = 'fdfmethod.lpnormlp1';
else    
    %#function fmethod.eqriplp
    designobj.equiripple = 'fmethod.eqriplp';
end

%#function fmethod.firlslp
designobj.firls = 'fmethod.firlslp';

if ~isempty(str)
    designobj = designobj.(str);
end

% [EOF]
