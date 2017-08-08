function designobj = getdesignobj(~, str, sigonlyflag)
%GETDESIGNOBJ   Get the designobj.

%   Copyright 1988-2013 The MathWorks, Inc.

if nargin < 2
    str = [];
    sigonlyflag = false;
elseif nargin < 3
    sigonlyflag = false;
end

%#function fmethod.eqripbs
%#function fmethod.firlsbs
designobj.equiripple = 'fmethod.eqripbs';
designobj.firls      = 'fmethod.firlsbs';

if isfdtbxinstalled && ~sigonlyflag
    %#function fdfmethod.lpnormbs1
    designobj.iirlpnorm  = 'fdfmethod.lpnormbs1';
    %#function fdfmethod.eqripbs
    designobj.equiripple = 'fdfmethod.eqripbs';
end

if ~isempty(str)
    designobj = designobj.(str);
end

% [EOF]
