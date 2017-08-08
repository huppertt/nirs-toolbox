function s = saveobj(this)
%SAVEOBJ   Save this object.

%   Author(s): J. Schickler
%   Copyright 2005 The MathWorks, Inc.

s.class      = class(this);
s.AllFDesign = get(this, 'AllFDesign');
s.Response   = get(this, 'Response');
s.ratechangefactors = getratechangefactors(this);

% [EOF]
