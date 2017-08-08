function cSpecCon = getconstructor(this, stype)
%GETCONSTRUCTOR   Return the constructor for the specification type.

%   Author(s): J. Schickler
%   Copyright 1999-2005 The MathWorks, Inc.

if nargin < 2
    stype = get(this, 'SpecificationType');
end

switch lower(stype)
    case 'fp,fst,ap,ast',
        %#function fspecs.lpmin
        cSpecCon = 'fspecs.lpmin';
    case 'n,f3db',
        %#function fspecs.lp3db
        cSpecCon = 'fspecs.lp3db';
    case 'nb,na,f3db'
        %#function fspecs.lpiir3db
        cSpecCon = 'fspecs.lpiir3db';
    case 'n,fc',
        %#function fspecs.lpcutoff
        cSpecCon = 'fspecs.lpcutoff';
    case 'n,f3db,ap',
        %#function fspecs.lpcutoffwap
        cSpecCon = 'fspecs.lpcutoffwap';
    case 'n,fp,f3db',
        %#function fspecs.lpcutoffwfp
        cSpecCon = 'fspecs.lpcutoffwfp';
    case 'n,f3db,ast',
        %#function fspecs.lpcutoffwas
        cSpecCon = 'fspecs.lpcutoffwas';
    case 'n,f3db,fst',
        %#function fspecs.lpcutoffwfs
        cSpecCon = 'fspecs.lpcutoffwfs';
    case 'n,fc,ap,ast',
        %#function fspecs.lpcutoffwatten
        cSpecCon = 'fspecs.lpcutoffwatten';
    case 'n,fp,ap',
        %#function fspecs.lppass
        cSpecCon = 'fspecs.lppass';
    case 'n,fst,ast',
        %#function fspecs.lpstop
        cSpecCon = 'fspecs.lpstop';
    case 'n,fp,fst,ast'
        %#function fspecs.lpstopfpass
        cSpecCon = 'fspecs.lpstopfpass';
    case 'n,fst,ap,ast'
        %#function fspecs.lpstopapass
        cSpecCon = 'fspecs.lpstopapass';
    case 'n,fp,ap,ast',
        %#function fspecs.lppassastop
        cSpecCon = 'fspecs.lppassastop';
    case 'n,f3db,ap,ast',
        %#function fspecs.lpcutoffwapas
        cSpecCon = 'fspecs.lpcutoffwapas';
    case 'n,fp,fst,ap',
        %#function fspecs.lppassfstop
        cSpecCon = 'fspecs.lppassfstop';
    case 'n,fp,fst'
        %#function fspecs.lpweight
        cSpecCon = 'fspecs.lpweight';
    case 'nb,na,fp,fst'
        %#function fspecs.lpiir
        cSpecCon = 'fspecs.lpiir';
    otherwise
        error(message('signal:fdesign:lowpass:getconstructor:internalError'));
end
                            

% [EOF]
