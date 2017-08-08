function s = saveobj(this)
%SAVEOBJ   Save this object.

%   Author(s): J. Schickler
%   Copyright 1988-2005 The MathWorks, Inc.

s = savemetadata(this);

s.class = class(this);
s.version = this.version;

s.PersistentMemory    = this.PersistentMemory;
s.NumSamplesProcessed = this.NumSamplesProcessed;

for indx = 1:nstages(this)
    s.Stage(indx) = this.Stage(indx);
end

% [EOF]
