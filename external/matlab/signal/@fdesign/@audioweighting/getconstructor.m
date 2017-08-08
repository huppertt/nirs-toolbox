function cSpecCon = getconstructor(this, stype)
%GETCONSTRUCTOR   Get the constructor.

%   Copyright 2009 The MathWorks, Inc.

if nargin < 2
    stype = get(this, 'SpecificationType');
end

switch lower(stype)
    case 'wt',
        %#function fspecs.audioweightingwt
        cSpecCon = 'fspecs.audioweightingwt';
    case 'wt,class',
        %#function fspecs.audioweightingwtclass
        cSpecCon = 'fspecs.audioweightingwtclass';        
    otherwise
        error(message('signal:fdesign:audioweighting:getconstructor:internalError'));
end


% [EOF]

