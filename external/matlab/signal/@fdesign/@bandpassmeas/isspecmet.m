function b = isspecmet(this, hfdesign)
%ISSPECMET   True if the object is specmet.

%   Author(s): J. Schickler
%   Copyright 2005 The MathWorks, Inc.

if nargin < 2
    hfdesign = get(this, 'Specification');
end

% Get the specifications from the object.
specs = measureinfo(hfdesign);

% Return true if the attenuations are greater than the specification and
% the ripples are smaller than the specification.
if ((isempty(specs.Astop1) || this.Astop1 >= specs.Astop1) && ...
    (isempty(specs.Apass)  || this.Apass  <= specs.Apass) && ...
    (isempty(specs.Astop2) || this.Astop2 >= specs.Astop2))
    b = true;
else
    b = false;
end


% [EOF]
