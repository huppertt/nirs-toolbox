function cSpecCon = getconstructor(this, stype)
%GETCONSTRUCTOR   Return the constructor for the specification type.

%   Author(s): J. Schickler
%   Copyright 1999-2005 The MathWorks, Inc.

if nargin < 2
    stype = get(this, 'Specification');
end

switch lower(stype)
    case 'tw,ast',
        if this.Band == 2
            %#function fspecs.hbmin
            cSpecCon = 'fspecs.hbmin';
        else
            %#function fspecs.nyqmin
            cSpecCon = 'fspecs.nyqmin';
        end
    case 'n,tw',
        if this.Band == 2
            %#function fspecs.hbordntw
            cSpecCon = 'fspecs.hbordntw';
        else
            %#function fspecs.nyqordntw
            cSpecCon = 'fspecs.nyqordntw';
        end
    case 'n',
        if this.Band == 2
            %#function fspecs.hbord
            cSpecCon = 'fspecs.hbord';
        else
            %#function fspecs.nyqord
            cSpecCon = 'fspecs.nyqord';
        end
    case 'n,ast',
        if this.Band == 2
            %#function fspecs.hbordastop
            cSpecCon = 'fspecs.hbordastop';
        else
            %#function fspecs.nyqordastop
            cSpecCon = 'fspecs.nyqordastop';
        end
    otherwise
        error(message('signal:fdesign:nyquist:getconstructor:internalError'));
end

% [EOF]
