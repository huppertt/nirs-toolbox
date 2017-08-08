function s = saveobj(this)
%SAVEOBJ   Save this object.

%   Author(s): J. Schickler
%   Copyright 1988-2004 The MathWorks, Inc.

s          = get(this);
s.Fs       = get(this, 'privFs');
s.class    = class(this);
s.Metadata = get(this, 'Metadata');
s.CenterDC = get(this, 'CenterDC');

s = setstructfields(s, thissaveobj(this));

% [EOF]
