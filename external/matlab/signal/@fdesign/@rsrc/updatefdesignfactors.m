function updatefdesignfactors(this)
%UPDATEFDESIGNFACTORS   Update the CurrentFDesign rate change factors.

%   Author(s): J. Schickler
%   Copyright 2005-2008 The MathWorks, Inc.

interpfactor = get(this, 'InterpolationFactor');
decimfactor  = get(this, 'DecimationFactor');

if isprop(this.CurrentFDesign, 'DecimationFactor')
    set(this.CurrentFDesign, 'DecimationFactor', decimfactor);
end
if isprop(this.CurrentFDesign, 'InterpolationFactor')
    set(this.CurrentFDesign, 'InterpolationFactor', interpfactor);
end

% [EOF]
