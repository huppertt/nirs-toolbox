function designobj = getdesignobj(~, str, sigonlyflag)
%GETDESIGNOBJ Get the designobj.

%   Copyright 1988-2013 The MathWorks, Inc.

if nargin < 2
    str = [];
    sigonlyflag = false;
elseif nargin < 3
    sigonlyflag = false;
end

%#function fmethod.eqripbp
%#function fmethod.firlsbp
designobj.equiripple = 'fmethod.eqripbp';
designobj.firls      = 'fmethod.firlsbp';

if isfdtbxinstalled && ~sigonlyflag
    %#function fdfmethod.lpnormbp1
    designobj.iirlpnorm  = 'fdfmethod.lpnormbp1';
    %#function fdfmethod.eqripbp
    designobj.equiripple = 'fdfmethod.eqripbp';
end

if ~isempty(str)
    designobj = designobj.(str);
end

% [EOF]
