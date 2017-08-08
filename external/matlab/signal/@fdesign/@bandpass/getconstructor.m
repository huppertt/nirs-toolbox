function cSpecCon = getconstructor(this, stype)
%GETCONSTRUCTOR Return the constructor for the specification type.

%   Copyright 1999-2011 The MathWorks, Inc.

if nargin < 2
    stype = get(this, 'SpecificationType');
end

switch lower(stype)
    case 'fst1,fp1,fp2,fst2,ast1,ap,ast2',
        %#function fspecs.bpmin
        cSpecCon = 'fspecs.bpmin';
    case 'n,f3db1,f3db2',
        %#function fspecs.bp3db
        cSpecCon = 'fspecs.bp3db';
    case 'n,fc1,fc2',
        %#function fspecs.bpcutoff
        cSpecCon = 'fspecs.bpcutoff';
    case 'n,fc1,fc2,ast1,ap,ast2'
        %#function fspecs.bpcutoffwatten
        cSpecCon = 'fspecs.bpcutoffwatten';
    case 'n,f3db1,f3db2,ap',
        %#function fspecs.bpcutoffwap
        cSpecCon = 'fspecs.bpcutoffwap';
    case 'n,f3db1,f3db2,bwp',
        %#function fspecs.bpcutoffwbwp
        cSpecCon = 'fspecs.bpcutoffwbwp';
    case 'n,f3db1,f3db2,ast',
        %#function fspecs.bpcutoffwas
        cSpecCon = 'fspecs.bpcutoffwas';
    case 'n,f3db1,f3db2,bwst',
        %#function fspecs.bpcutoffwbws
        cSpecCon = 'fspecs.bpcutoffwbws';
    case 'n,fp1,fp2,ap',
        %#function fspecs.bppass
        cSpecCon = 'fspecs.bppass';
    case 'n,fst1,fst2,ast',
        %#function fspecs.bpstop
        cSpecCon = 'fspecs.bpstop';
    case 'n,fp1,fp2,ast1,ap,ast2',
        %#function fspecs.bppassastop
        cSpecCon = 'fspecs.bppassastop';
    case 'n,f3db1,f3db2,ast1,ap,ast2',
        %#function fspecs.bpcutoffwapas
        cSpecCon = 'fspecs.bpcutoffwapas';
    case 'n,fst1,fp1,fp2,fst2,ap',
        %#function fspecs.bppassfstop
        cSpecCon = 'fspecs.bppassfstop';
    case 'n,fst1,fp1,fp2,fst2'
        %#function fspecs.bpweight
        cSpecCon = 'fspecs.bpweight';
    case 'n,fst1,fp1,fp2,fst2,c'
        %#function fspecs.bpweightconstrained
        cSpecCon = 'fspecs.bpweightconstrained';        
    case 'nb,na,fst1,fp1,fp2,fst2'
        %#function fspecs.bpiir
        cSpecCon = 'fspecs.bpiir';
    otherwise
        error(message('signal:fdesign:bandpass:getconstructor:internalError'));
end

% [EOF]
