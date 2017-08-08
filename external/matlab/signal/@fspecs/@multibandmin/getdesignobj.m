function designobj = getdesignobj(~, str)
%GETDESIGNOBJ Get the design object.

%   Copyright 2011 The MathWorks, Inc.

%#function fdfmethod.eqripmultibandmin
designobj.equiripple = 'fdfmethod.eqripmultibandmin';

if nargin > 1
    designobj = designobj.(str);
end

% [EOF]
