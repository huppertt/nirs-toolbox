function c = getconstructor(this, stype)
%GETCONSTRUCTOR   Get the constructor.

%   Author(s): J. Schickler
%   Copyright 2005 The MathWorks, Inc.

if nargin < 2
    stype = get(this, 'SpecificationType');
end

switch lower(stype)
    case 'fp,fst,ap,ast'
        %#function fspecs.lpisincmin
        c = 'fspecs.lpisincmin';
    case 'n,fc,ap,ast'
        %#function fspecs.lpisinccutoffwatten
        c = 'fspecs.lpisinccutoffwatten';
    case 'n,fp,ap,ast'
        %#function fspecs.lpisincpassastop
        c = 'fspecs.lpisincpassastop';
    case 'n,fst,ap,ast'
        %#function fspecs.lpisincstopapass
        c = 'fspecs.lpisincstopapass';
    case 'n,fp,fst'
        %#function fspecs.lpisinc
        c = 'fspecs.lpisinc';
end

% [EOF]
