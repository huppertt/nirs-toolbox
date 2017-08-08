function designobj = getdesignobj(~, str)
%GETDESIGNOBJ Get the design object.

%   Copyright 2011 The MathWorks, Inc.

%#function fdfmethod.eqriphpminisinc
designobj.equiripple = 'fdfmethod.eqriphpminisinc';

if nargin > 1
  designobj = designobj.(str);
end

% [EOF]
