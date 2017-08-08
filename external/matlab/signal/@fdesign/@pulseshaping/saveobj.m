function s = saveobj(this)
%SAVEOBJ   Save this object.

%   Copyright 2008 The MathWorks, Inc.

s.class         = class(this);
s.Response      = get(this, 'Response');
s.PulseShape    = get(this, 'PulseShape');
s.PulseShapeObj = saveobj(this.PulseShapeObj);
s.version       = '9a';

% [EOF]
