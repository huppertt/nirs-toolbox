function cSpecCon = getconstructor(this, stype)
%GETCONSTRUCTOR   Get the constructor.

%   Author(s): V. Pellissier
%   Copyright 2005 The MathWorks, Inc.

if nargin < 2
    stype = get(this, 'SpecificationType');
end

switch lower(stype)
    case 'n,f0',
        %#function fspecs.octavewordncenterfreq
        cSpecCon = 'fspecs.octavewordncenterfreq';
    otherwise
        error(message('signal:fdesign:octave:getconstructor:internalError'));
end


% [EOF]

