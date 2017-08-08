function ratechangefactors = getratechangefactors(this)
%GETRATECHANGEFACTORS   Get the ratechangefactors.

%   Author(s): J. Schickler
%   Copyright 2005 The MathWorks, Inc.

ratechangefactors = [get(this, 'InterpolationFactor') 1];

% [EOF]
