function cSpecCon = getconstructor(this, stype)
%GETCONSTRUCTOR Get the constructor.

%   Copyright 2005-2011 The MathWorks, Inc.

if nargin < 2
    stype = get(this, 'SpecificationType');
end

switch lower(stype)
    case 'n,f,a'
        %#function fspecs.sbarbmag
        cSpecCon = 'fspecs.sbarbmag';
    case 'f,a,r'
         %#function fspecs.sbarbmagmin
         cSpecCon = 'fspecs.sbarbmagmin';
    case 'nb,na,f,a'
        %#function fspecs.sbarbmagiir
        cSpecCon = 'fspecs.sbarbmagiir';
    case 'n,b,f,a'
        %#function fspecs.multiband
        cSpecCon = 'fspecs.multiband';
    case 'n,b,f,a,c'
        %#function fspecs.multibandconstrained
        cSpecCon = 'fspecs.multibandconstrained';        
    case 'b,f,a,r'
        %#function = fspecs.multibandmin
        cSpecCon = 'fspecs.multibandmin';
    case 'nb,na,b,f,a'
        %#function fspecs.multibandiir
        cSpecCon = 'fspecs.multibandiir';
    otherwise
        error(message('signal:fdesign:arbmag:getconstructor:internalError'));
end


% [EOF]
