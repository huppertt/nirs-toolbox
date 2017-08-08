function designobj = getdesignobj(~, str, sigonlyflag)
%GETDESIGNOBJ   Get the designobj.

%   Copyright 2008-2013 The MathWorks, Inc.

if nargin < 2
    str = [];
    sigonlyflag = false;
elseif nargin < 3
    sigonlyflag = false;
end

%#function fmethod.firlsdifford
designobj.firls = 'fmethod.firlsdifford';

if isfdtbxinstalled && ~sigonlyflag
   %#function fdfmethod.eqripdifford
    designobj.equiripple = 'fdfmethod.eqripdifford';
else
   %#function fmethod.eqripdifford
    designobj.equiripple = 'fmethod.eqripdifford';
end

if ~isempty(str)
    designobj = designobj.(str);
end

% [EOF]
