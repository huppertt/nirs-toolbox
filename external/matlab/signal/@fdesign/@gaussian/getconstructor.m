function cSpecCon = getconstructor(this)
%GETCONSTRUCTOR Return the constructor for the specification type.

%   Copyright 2008 The MathWorks, Inc.


if nargin < 2
    stype = get(this, 'SpecificationType');
end

switch lower(stype)
    case 'nsym,bt',
        %#function fspecs.psgaussnsym
        cSpecCon = 'fspecs.psgaussnsym';
    otherwise
        error(message('signal:fdesign:gaussian:getconstructor:internalError'));
end

% [EOF]
