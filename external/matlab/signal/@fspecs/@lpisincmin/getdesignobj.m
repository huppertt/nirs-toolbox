function designobj = getdesignobj(this, str)
%GETDESIGNOBJ Get the design object.

%   Copyright 2005-2011 The MathWorks, Inc.

if this.privInvSincParamsTunableFlag
  %#function fdfmethod.eqriplpminisincwparams
  designobj.equiripple = 'fdfmethod.eqriplpminisincwparams';
else
  %#function fdfmethod.eqriplpminisinc
  designobj.equiripple = 'fdfmethod.eqriplpminisinc';
end
if nargin > 1
  designobj = designobj.(str);
end

% [EOF]
