function designobj = getdesignobj(~, str, sigonlyflag)
%GETDESIGNOBJ Get the designobj.

%   Copyright 2008-2013 The MathWorks, Inc.

if nargin < 2
    str = [];
    sigonlyflag = false;
elseif nargin < 3
    sigonlyflag = false;
end

%#function fmethod.firlsdiffordmb
designobj.firls = 'fmethod.firlsdiffordmb';

if isfdtbxinstalled && ~sigonlyflag
   %#function fdfmethod.eqripdiffordmb
    designobj.equiripple = 'fdfmethod.eqripdiffordmb';
else
   %#function fmethod.eqripdiffordmb
    designobj.equiripple = 'fmethod.eqripdiffordmb';
end

if ~isempty(str)
    designobj = designobj.(str);
end

% [EOF]
