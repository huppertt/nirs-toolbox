function designobj = getdesignobj(this, str)
%GETDESIGNOBJ Get the design object.

%   Copyright 2005-2011 The MathWorks, Inc.

if this.privInvSincParamsTunableFlag
  %#function fdfmethod.eqriplpapassisincwparams
  designobj.equiripple = 'fdfmethod.eqriplpapassisincwparams';
else
  %#function fdfmethod.eqriplpapassisinc
  designobj.equiripple = 'fdfmethod.eqriplpapassisinc';
end

if nargin > 1
  designobj = designobj.(str);
end

% [EOF]
