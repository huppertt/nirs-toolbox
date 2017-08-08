function cSpecCon = getconstructor(this, stype)
%GETCONSTRUCTOR Return the constructor for the specification type.

%   Copyright 1999-2011 The MathWorks, Inc.

if nargin < 2
    stype = get(this, 'SpecificationType');
end

switch lower(stype)
    case 'fp1,fst1,fst2,fp2,ap1,ast,ap2',
        %#function fspecs.bsmin
        cSpecCon = 'fspecs.bsmin';
    case 'n,f3db1,f3db2',
        %#function fspecs.bs3db
        cSpecCon = 'fspecs.bs3db';
    case 'n,fc1,fc2',
        %#function fspecs.bscutoff
        cSpecCon = 'fspecs.bscutoff';
    case 'n,fc1,fc2,ap1,ast,ap2',
        %#function fspecs.bscutoffwatten
        cSpecCon = 'fspecs.bscutoffwatten';
    case 'n,f3db1,f3db2,ap',
        %#function fspecs.bscutoffwap
        cSpecCon = 'fspecs.bscutoffwap';
    case 'n,f3db1,f3db2,bwp',
        %#function fspecs.bscutoffwbwp
        cSpecCon = 'fspecs.bscutoffwbwp';
    case 'n,f3db1,f3db2,ast',
        %#function fspecs.bscutoffwas
        cSpecCon = 'fspecs.bscutoffwas';
    case 'n,f3db1,f3db2,bwst',
        %#function fspecs.bscutoffwbws
        cSpecCon = 'fspecs.bscutoffwbws';
    case 'n,fp1,fp2,ap',
        %#function fspecs.bspass
        cSpecCon = 'fspecs.bspass';
    case 'n,fst1,fst2,ast',
        %#function fspecs.bsstop
        cSpecCon = 'fspecs.bsstop';
    case 'n,fp1,fp2,ap,ast',
        %#function fspecs.bspassastop
        cSpecCon = 'fspecs.bspassastop';
    case 'n,f3db1,f3db2,ap,ast',
        %#function fspecs.bscutoffwapas
        cSpecCon = 'fspecs.bscutoffwapas';
    case 'n,fp1,fst1,fst2,fp2,ap',
        %#function fspecs.bspassfstop
        cSpecCon = 'fspecs.bspassfstop';
    case 'n,fp1,fst1,fst2,fp2'
        %#function fspecs.bsweight
        cSpecCon = 'fspecs.bsweight';
    case 'n,fp1,fst1,fst2,fp2,c'
        %#function fspecs.bsweightconstrained
        cSpecCon = 'fspecs.bsweightconstrained';                
    case 'nb,na,fp1,fst1,fst2,fp2'
        %#function fspecs.bsiir
        cSpecCon = 'fspecs.bsiir';
        
    otherwise
        error(message('signal:fdesign:bandstop:getconstructor:internalError'));
end

% [EOF]
