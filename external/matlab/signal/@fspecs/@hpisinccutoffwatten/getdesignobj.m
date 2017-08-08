function designobj = getdesignobj(~, str)
%GETDESIGNOBJ Get the design object.

%   Copyright 2011 The MathWorks, Inc.

%#function fdfmethod.eqriphpcutoffisinc
designobj.equiripple = 'fdfmethod.eqriphpcutoffisinc';

if nargin > 1
    designobj = designobj.(str);
end

% [EOF]
