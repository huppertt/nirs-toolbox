function designobj = getdesignobj(this, str)
%GETDESIGNOBJ   Get the design object.

%   Author(s): J. Schickler
%   Copyright 1988-2005 The MathWorks, Inc.

%#function fmethod.elliplpastop
designobj.ellip = 'fmethod.elliplpastop';

if nargin > 1
    designobj = designobj.(str);
end

% [EOF]
