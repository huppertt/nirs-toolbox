function designobj = getdesignobj(~, str, sigonlyflag)
%GETDESIGNOBJ   Get the designobj.

%   Copyright 1988-2013 The MathWorks, Inc.

if nargin < 2
    str = [];
    sigonlyflag = false;
elseif nargin < 3
    sigonlyflag = false;
end

%#function fmethod.elliphpastop
designobj.ellip      = 'fmethod.elliphpastop';

if isfdtbxinstalled && ~sigonlyflag
    %#function fdfmethod.eqriphpastop
    designobj.equiripple = 'fdfmethod.eqriphpastop';
end

if ~isempty(str)
    designobj = designobj.(str);
end

% [EOF]
