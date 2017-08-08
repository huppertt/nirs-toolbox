function cSpecCon = getconstructor(this, stype)
%GETCONSTRUCTOR Return the constructor for the specification type.

%   Copyright 1999-2011 The MathWorks, Inc.

if nargin < 2
    stype = get(this, 'SpecificationType');
end

switch lower(stype)
    case 'fst,fp,ast,ap',
        %#function fspecs.hpmin
        cSpecCon = 'fspecs.hpmin';
    case 'n,f3db',
        %#function fspecs.hp3db
        cSpecCon = 'fspecs.hp3db';
    case 'nb,na,f3db',
        %#function fspecs.hpiir3db
        cSpecCon = 'fspecs.hpiir3db';        
    case 'n,fc',                
        %#function fspecs.hpcutoff
        cSpecCon = 'fspecs.hpcutoff';
    case 'n,f3db,ap',
        %#function fspecs.hpcutoffwap
        cSpecCon = 'fspecs.hpcutoffwap';
    case 'n,f3db,fp',
        %#function fspecs.hpcutoffwfp
        cSpecCon = 'fspecs.hpcutoffwfp';
    case 'n,f3db,ast',
        %#function fspecs.hpcutoffwas
        cSpecCon = 'fspecs.hpcutoffwas';
    case 'n,fc,ast,ap',
        %#function fspecs.hpcutoffwatten
        cSpecCon = 'fspecs.hpcutoffwatten';
    case 'n,fst,f3db',
        %#function fspecs.hpcutoffwfs
        cSpecCon = 'fspecs.hpcutoffwfs';
    case 'n,fp,ap',
        %#function fspecs.hppass
        cSpecCon = 'fspecs.hppass';
    case 'n,fp,ast,ap',
        %#function fspecs.hppassastop
        cSpecCon = 'fspecs.hppassastop';
    case 'n,f3db,ast,ap',
        %#function fspecs.hpcutoffwapas
        cSpecCon = 'fspecs.hpcutoffwapas';
    case 'n,fst,ast',
        %#function fspecs.hpstop
        cSpecCon = 'fspecs.hpstop';
    case 'n,fst,ast,ap',
        %#function fspecs.hpstopapass
        cSpecCon = 'fspecs.hpstopapass';
    case 'n,fst,fp,ast',
        %#function fspecs.hpstopfpass
        cSpecCon = 'fspecs.hpstopfpass';
    case 'n,fst,fp,ap',
        %#function fspecs.hppassfstop
        cSpecCon = 'fspecs.hppassfstop';
    case 'n,fst,fp'
        %#function fspecs.hpweight
        cSpecCon = 'fspecs.hpweight';
    case 'nb,na,fst,fp'
        %#function fspecs.hpiir
        cSpecCon = 'fspecs.hpiir';
    otherwise
        error(message('signal:fdesign:highpass:getconstructor:internalError'));
end

% [EOF]
