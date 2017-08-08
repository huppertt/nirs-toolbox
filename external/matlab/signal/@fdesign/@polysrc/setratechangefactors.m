function setratechangefactors(this, ratechangefactors)
%SETRATECHANGEFACTORS   Set the ratechangefactors.

%   Author(s): J. Schickler
%   Copyright 2005 The MathWorks, Inc.

set(this, 'InterpolationFactor', ratechangefactors(1), ...
    'DecimationFactor', ratechangefactors(2));

% [EOF]
