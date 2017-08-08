function cSpecCon = getconstructor(this, stype)
%GETCONSTRUCTOR   Return the constructor for the specification type.

%   Author(s): R. Losada
%   Copyright 2006 The MathWorks, Inc.

if nargin < 2
    stype = get(this, 'SpecificationType');
end

switch lower(stype)
    case 'n,f0,q',
        %#function fspecs.peakq
        cSpecCon = 'fspecs.peakq'; 
    case 'n,f0,q,ap',
        %#function fspecs.peakqap
        cSpecCon = 'fspecs.peakqap'; 
    case 'n,f0,q,ast',
        %#function fspecs.peakqast
        cSpecCon = 'fspecs.peakqast';         
    case 'n,f0,q,ap,ast',
        %#function fspecs.peakqapast
        cSpecCon = 'fspecs.peakqapast';                        
    case 'n,f0,bw',
        %#function fspecs.peakbw
        cSpecCon = 'fspecs.peakbw';                        
    case 'n,f0,bw,ap',
        %#function fspecs.peakbwap
        cSpecCon = 'fspecs.peakbwap';                        
    case 'n,f0,bw,ast',
        %#function fspecs.peakbwast
        cSpecCon = 'fspecs.peakbwast';                        
    case 'n,f0,bw,ap,ast',
        %#function fspecs.peakbwapast
        cSpecCon = 'fspecs.peakbwapast';                
    otherwise
        error(message('signal:fdesign:peak:getconstructor:internalError'));
end

% [EOF]
