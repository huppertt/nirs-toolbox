function ratechangefactors = getratechangefactors(this)
%GETRATECHANGEFACTORS   Get the ratechangefactors.

%   Author(s): J. Schickler
%   Copyright 2005 The MathWorks, Inc.

ratechangefactors = [1 get(this, 'DecimationFactor')];

% [EOF]
