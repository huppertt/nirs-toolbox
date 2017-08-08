function disp(this)
%DISP   Display this object.

%   Author(s): J. Schickler
%   Copyright 1988-2004 The MathWorks, Inc.

switch lower(this.View)
    case {'complete', 'individual'}
        props = {'View'};
    case 'cumulative'
        props = {'View', 'SecondaryScaling'};
    case 'userdefined'
        props = {'View', 'UserDefinedSections'};
end

siguddutils('dispstr', this, props);

% [EOF]
