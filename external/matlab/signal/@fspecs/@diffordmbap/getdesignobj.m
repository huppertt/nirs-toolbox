function designobj = getdesignobj(~, str)
%GETDESIGNOBJ Get the designobj.

%   Copyright 2011 The MathWorks, Inc.

%#function fdfmethod.eqripdiffordmbap
designobj.equiripple = 'fdfmethod.eqripdiffordmbap';

if nargin > 1
  designobj = designobj.(str);
end

% [EOF]
