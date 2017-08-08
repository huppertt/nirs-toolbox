function designobj = getdesignobj(this, str)
%GETDESIGNOBJ Get the design object.

%   Copyright 2005-2011 The MathWorks, Inc.

if this.privInvSincParamsTunableFlag
  %#function fdfmethod.eqriplpastopisincwparams
  designobj.equiripple = 'fdfmethod.eqriplpastopisincwparams';
else
  %#function fdfmethod.eqriplpastopisinc
  designobj.equiripple = 'fdfmethod.eqriplpastopisinc';
end
if nargin > 1
  designobj = designobj.(str);
end

% [EOF]
