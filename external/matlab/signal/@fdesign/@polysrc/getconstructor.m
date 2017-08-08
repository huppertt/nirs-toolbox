function cSpecCon = getconstructor(this, stype)
%GETCONSTRUCTOR   Get the constructor.

%   Copyright 2007 The MathWorks, Inc.

if nargin < 2
    stype = get(this, 'SpecificationType');
end

switch lower(stype)
    case 'np',
        %#function fspecs.fdsrcword
        cSpecCon = 'fspecs.fdsrcword';
    otherwise
        error(message('signal:fdesign:polysrc:getconstructor:internalError'));
end


% [EOF]

