function designobj = getdesignobj(~, str, sigonlyflag)
%GETDESIGNOBJ   Get the designobj.

%   Copyright 1988-2013 The MathWorks, Inc.

if nargin < 2
    str = [];
    sigonlyflag = false;
elseif nargin < 3
    sigonlyflag = false;
end

%#function fmethod.eqriphp
designobj.equiripple = 'fmethod.eqriphp';
%#function fmethod.firlshp
designobj.firls      = 'fmethod.firlshp';

if isfdtbxinstalled && ~sigonlyflag
    %#function fdfmethod.lpnormhp1
    designobj.iirlpnorm = 'fdfmethod.lpnormhp1';
    %#function fdfmethod.eqriphp
    designobj.equiripple = 'fdfmethod.eqriphp';
end

if ~isempty(str)
    designobj = designobj.(str);
end

% [EOF]
