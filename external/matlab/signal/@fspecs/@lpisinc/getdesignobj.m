function designobj = getdesignobj(this, str)
%GETDESIGNOBJ Get the design object.

%   Copyright 2005-2011 The MathWorks, Inc.

if this.privInvSincParamsTunableFlag
  %#function fdfmethod.eqriplpisincwparams
  designobj.equiripple = 'fdfmethod.eqriplpisincwparams';
else
  %#function fdfmethod.eqriplpisinc
  designobj.equiripple = 'fdfmethod.eqriplpisinc';
end

if nargin > 1
  designobj = designobj.(str);
end

% [EOF]
