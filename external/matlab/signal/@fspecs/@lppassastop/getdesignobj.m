function designobj = getdesignobj(~, str, sigonlyflag)
%GETDESIGNOBJ   Get the designobj.

%   Copyright 1988-2013 The MathWorks, Inc.

if nargin < 2
    str = [];
    sigonlyflag = false;
elseif nargin < 3
    sigonlyflag = false;
end

%#function fmethod.elliplpastop
designobj.ellip      = 'fmethod.elliplpastop';

if isfdtbxinstalled && ~sigonlyflag
    %#function fdfmethod.eqriplpastop
    designobj.equiripple = 'fdfmethod.eqriplpastop';
end

if ~isempty(str)
    designobj = designobj.(str);
end

% [EOF]
