function cSpecCon = getconstructor(this, stype)
%GETCONSTRUCTOR Return the constructor for the specification type.

%   Copyright 2005-2011 The MathWorks, Inc.

if nargin < 2
    stype = get(this, 'SpecificationType');
end

switch lower(stype)
    case 'n',
        %#function fspecs.difford
        cSpecCon = 'fspecs.difford';      % Specify order, single-band
    case 'n,fp,fst',
        %#function fspecs.diffordmb
        cSpecCon = 'fspecs.diffordmb';    % Specify order, multi-band
    case 'n,fp,fst,ap',
        %#function fspecs.diffordmbap
        cSpecCon = 'fspecs.diffordmbap';  % Specify order and Apass, multi-band
    case 'n,fp,fst,ast',
        %#function fspecs.diffordmbast
        cSpecCon = 'fspecs.diffordmbast';  % Specify order and Astop, multi-band        
    case 'ap',
        %#function fspecs.diffmin
        cSpecCon = 'fspecs.diffmin';      % Minimum-order, single-band
    case 'fp,fst,ap,ast',
        %#function fspecs.diffminmb
        cSpecCon = 'fspecs.diffminmb';    % Minimum-order, multi-band
    otherwise
        error(message('signal:fdesign:differentiator:getconstructor:internalError'));
end

% [EOF]
