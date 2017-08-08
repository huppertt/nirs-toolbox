function updatefdesignfactors(this)
%UPDATEFDESIGNFACTORS   Update the current FDesign rate change factors.

%   Author(s): J. Schickler
%   Copyright 2005-2008 The MathWorks, Inc.

interpolationfactor = get(this, 'InterpolationFactor');

if isprop(this.CurrentFDesign, 'InterpolationFactor')
    set(this.CurrentFDesign, 'InterpolationFactor', interpolationfactor);
end

% [EOF]
