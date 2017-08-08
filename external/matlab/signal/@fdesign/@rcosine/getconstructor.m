function cSpecCon = getconstructor(this)
%GETCONSTRUCTOR Return the constructor for the specification type.

%   Copyright 2008 The MathWorks, Inc.


if nargin < 2
    stype = get(this, 'SpecificationType');
end

switch lower(stype)
    case 'nsym,beta',
        %#function fspecs.psrcosnsym
        cSpecCon = 'fspecs.psrcosnsym';
    case 'n,beta',
        %#function fspecs.psrcosord
        cSpecCon = 'fspecs.psrcosord';
    case 'ast,beta',
        %#function fspecs.psrcosmin
        cSpecCon = 'fspecs.psrcosmin';
    otherwise
        error(message('signal:fdesign:rcosine:getconstructor:internalError'));
end

% [EOF]
