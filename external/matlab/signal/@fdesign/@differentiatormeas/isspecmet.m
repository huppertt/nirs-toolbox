function b = isspecmet(this, hfdesign)
%ISSPECMET   True if the object is specmet.

%   Author(s): J. Schickler
%   Copyright 2005 The MathWorks, Inc.

if nargin < 2
    hfdesign = get(this, 'Specification');
end

% Get the specifications from the FDesign object.
specs = measureinfo(hfdesign);

% Return true if the measured Apass is less than or equal to the specificed
% and if the measured Astop is greater than or eual to the specified.
if  (isempty(specs.Apass) || this.Apass <= specs.Apass) && ...
    (isempty(specs.Astop) || this.Astop >= specs.Astop)
    b = true;
else
    b = false;
end

% [EOF]
