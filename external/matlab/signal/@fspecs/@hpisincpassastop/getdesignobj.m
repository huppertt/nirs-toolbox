function designobj = getdesignobj(~, str)
%GETDESIGNOBJ Get the design object.

%   Copyright 2011 The MathWorks, Inc.

%#function fdfmethod.eqriphpastopisinc
designobj.equiripple = 'fdfmethod.eqriphpastopisinc';

if nargin > 1
  designobj = designobj.(str);
end

% [EOF]
