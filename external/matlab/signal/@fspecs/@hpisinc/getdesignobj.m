function designobj = getdesignobj(~, str)
%GETDESIGNOBJ Get the design object.

%   Copyright 2011 The MathWorks, Inc.

%#function fdfmethod.eqriphpisinc
designobj.equiripple = 'fdfmethod.eqriphpisinc';

if nargin > 1
  designobj = designobj.(str);
end

% [EOF]
