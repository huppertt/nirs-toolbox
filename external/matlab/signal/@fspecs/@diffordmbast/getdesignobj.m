function designobj = getdesignobj(~, str)
%GETDESIGNOBJ Get the designobj.

%   Copyright 2011 The MathWorks, Inc.

%#function fdfmethod.eqripdiffordmbast
designobj.equiripple = 'fdfmethod.eqripdiffordmbast';

if nargin > 1
  designobj = designobj.(str);
end

% [EOF]
