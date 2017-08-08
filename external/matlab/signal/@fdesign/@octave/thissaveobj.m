function s = thissaveobj(this)
%THISSAVEOBJ   Save this object.

%   Author(s): J. Schickler
%   Copyright 1999-2005 The MathWorks, Inc.

s.BandsPerOctave = get(this, 'BandsPerOctave');
s.Mask = get(this, 'Mask');

% [EOF]
