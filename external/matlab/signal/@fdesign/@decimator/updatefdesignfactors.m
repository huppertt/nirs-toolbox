function updatefdesignfactors(this)
%UPDATEFDESIGNFACTORS   Update the CurrentFDesign ratechange factors.

%   Author(s): J. Schickler
%   Copyright 2005-2008 The MathWorks, Inc.

decimationfactor = get(this, 'DecimationFactor');

if isprop(this.CurrentFDesign, 'DecimationFactor')
    set(this.CurrentFDesign, 'DecimationFactor', decimationfactor);
end

% [EOF]
