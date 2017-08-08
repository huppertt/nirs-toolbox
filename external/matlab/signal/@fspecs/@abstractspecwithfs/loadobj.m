function this = loadobj(this, s)
%LOADOBJ   Load this object.

%   Author(s): J. Schickler
%   Copyright 1988-2005 The MathWorks, Inc.

if nargin < 2
    s    = this;
    this = feval(s.class);
end

% We can't set the Fs directly.
this.privFs = s.Fs;

% We can't set the NormalizedFrequency directly.
this.privNormalizedFreq = s.NormalizedFrequency;

s = rmfield(s, {'Fs', 'class','NormalizedFrequency'});

set(this, s);

% [EOF]
