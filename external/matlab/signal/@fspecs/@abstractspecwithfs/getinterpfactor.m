function interpfactor = getinterpfactor(this)
%GETINTERPFACTOR   Get the interpfactor.

%   Author(s): J. Schickler
%   Copyright 2005 The MathWorks, Inc.

if isprop(this, 'InterpolationFactor')
    interpfactor = get(this, 'InterpolationFactor');
elseif isprop(this, 'Band')
    interpfactor = get(this, 'Band');
else
    % We should probably error here instead, but this will help the
    % HALFBAND case, which does not have a property for the band because it
    % is always 2.
    interpfactor = 2;
end

% [EOF]
