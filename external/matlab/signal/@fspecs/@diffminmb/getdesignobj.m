function designobj = getdesignobj(this, str)
%GETDESIGNOBJ   Get the designobj.

%   Author(s): P. Costa
%   Copyright 2005 The MathWorks, Inc.

%#function fdfmethod.eqripdiffminmb
designobj.equiripple = 'fdfmethod.eqripdiffminmb';

if nargin > 1
    designobj = designobj.(str);
end

% [EOF]
