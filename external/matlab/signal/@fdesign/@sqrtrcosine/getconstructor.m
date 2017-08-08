function cSpecCon = getconstructor(this)
%GETCONSTRUCTOR Return the constructor for the specification type.

%   Copyright 2008 The MathWorks, Inc.


if nargin < 2
    stype = get(this, 'SpecificationType');
end

switch lower(stype)
    case 'nsym,beta',
        %#function fspecs.pssqrtrcosnsym
        cSpecCon = 'fspecs.pssqrtrcosnsym';
    case 'n,beta',
        %#function fspecs.pssqrtrcosord
        cSpecCon = 'fspecs.pssqrtrcosord';
    case 'ast,beta',
        %#function fspecs.pssqrtrcosmin
        cSpecCon = 'fspecs.pssqrtrcosmin';
    otherwise
        error(message('signal:fdesign:sqrtrcosine:getconstructor:internalError'));
end

% [EOF]
