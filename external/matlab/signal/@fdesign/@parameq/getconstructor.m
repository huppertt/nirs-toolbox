function cSpecCon = getconstructor(this, stype)
%GETCONSTRUCTOR   Return the constructor for the specification type.

%   Copyright 2008 The MathWorks, Inc.

if nargin < 2
    stype = get(this, 'SpecificationType');
end

switch lower(stype)
    case 'f0,bw,bwp,gref,g0,gbw,gp',
        %#function fspecs.parameqbwgbwap
        cSpecCon = 'fspecs.parameqbwgbwap'; 
    case 'f0,bw,bwst,gref,g0,gbw,gst',
        %#function fspecs.parameqbwgbwast
        cSpecCon = 'fspecs.parameqbwgbwast'; 
    case 'f0,bw,bwp,gref,g0,gbw,gp,gst',
        %#function fspecs.parameqbwgbwapast
        cSpecCon = 'fspecs.parameqbwgbwapast';         
    case 'n,f0,bw,gref,g0,gbw',
        %#function fspecs.parameq
        cSpecCon = 'fspecs.parameq';                        
    case 'n,f0,bw,gref,g0,gbw,gp',
        %#function fspecs.parameqap
        cSpecCon = 'fspecs.parameqap';                        
    case 'n,f0,bw,gref,g0,gbw,gst',
        %#function fspecs.parameqast
        cSpecCon = 'fspecs.parameqast';                        
    case 'n,f0,bw,gref,g0,gbw,gp,gst',
        %#function fspecs.parameqapast
        cSpecCon = 'fspecs.parameqapast';                        
    case 'n,flow,fhigh,gref,g0,gbw',
        %#function fspecs.parameqflfh
        cSpecCon = 'fspecs.parameqflfh';                
    case 'n,flow,fhigh,gref,g0,gbw,gp',
        %#function fspecs.parameqflfhap
        cSpecCon = 'fspecs.parameqflfhap';        
    case 'n,flow,fhigh,gref,g0,gbw,gst',
        %#function fspecs.parameqflfhast
        cSpecCon = 'fspecs.parameqflfhast';
    case 'n,flow,fhigh,gref,g0,gbw,gp,gst',
        %#function fspecs.parameqflfhapast
        cSpecCon = 'fspecs.parameqflfhapast';
    case 'n,f0,qa,gref,g0',
        %#function fspecs.parameqaudioqa
        cSpecCon = 'fspecs.parameqaudioqa';
     case 'n,f0,fc,qa,g0',
        %#function fspecs.parameqaudioshelfqa
        cSpecCon = 'fspecs.parameqaudioshelfqa';
     case 'n,f0,fc,s,g0',
        %#function fspecs.parameqaudioshelfs
        cSpecCon = 'fspecs.parameqaudioshelfs';       
    otherwise
        error(message('signal:fdesign:parameq:getconstructor:internalError'));
end

% [EOF]
