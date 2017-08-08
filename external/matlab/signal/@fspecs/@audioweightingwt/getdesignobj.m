function designobj = getdesignobj(this, str)
%GETDESIGNOBJ   Get the design object.

%   Copyright 2009 The MathWorks, Inc.

switch lower(this.WeightingType)
    case 'cmessage'
        %#function fdfmethod.freqsampaudioweightcmessage
        designobj.freqsamp = 'fdfmethod.freqsampaudioweightcmessage';
        %#function fdfmethod.eqripaudioweightcmessage
        designobj.equiripple = 'fdfmethod.eqripaudioweightcmessage';
        %#function fdfmethod.bell41099audioweight
        designobj.bell41009 = 'fdfmethod.bell41099audioweight';
    case 'itut041'
        %#function fdfmethod.freqsampaudioweightitut041
        designobj.freqsamp = 'fdfmethod.freqsampaudioweightitut041';
        %#function fdfmethod.eqripaudioweightitut041
        designobj.equiripple = 'fdfmethod.eqripaudioweightitut041';
    case 'itur4684'
        %#function fdfmethod.freqsampaudioweightitur4684
        designobj.freqsamp = 'fdfmethod.freqsampaudioweightitur4684';
        %#function fdfmethod.eqripaudioweightitur4684
        designobj.equiripple = 'fdfmethod.eqripaudioweightitur4684';
        %#function fdfmethod.lpnormaudioweightitur4684
        designobj.iirlpnorm  = 'fdfmethod.lpnormaudioweightitur4684';
end

if nargin > 1
    designobj = designobj.(str);
end

% [EOF]
