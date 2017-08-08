function cSpecCon = getconstructor(this, stype)
%GETCONSTRUCTOR   Return the constructor for the specification type.

%   Author(s): P. Costa
%   Copyright 2005 The MathWorks, Inc.

if nargin < 2
    stype = get(this, 'SpecificationType');
end

switch lower(stype)
    case 'n,tw',
        %#function fspecs.hilbord
        cSpecCon = 'fspecs.hilbord';
    case 'tw,ap',
        %#function fspecs.hilbmin
        cSpecCon = 'fspecs.hilbmin';
    otherwise
        error(message('signal:fdesign:hilbert:getconstructor:internalError'));
end

% [EOF]
