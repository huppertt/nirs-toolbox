function constructor = getconstructor(this, dtype)
%GETCONSTRUCTOR Get the constructor.

%   Copyright 2005-2011 The MathWorks, Inc.

if nargin < 2
    dtype = get(this, 'Response');
end

switch lower(dtype)
    case 'inverse-sinc lowpass'
        %#function fdesign.isinclp
        constructor = 'isinclp';
    case 'inverse-sinc highpass'
        %#function fdesign.isinchp
        constructor = 'isinchp';        
    case 'cic compensator'
        %#function fdesign.ciccomp
        constructor = 'ciccomp';
    case 'cic'
        constructor = getcicconstructor(this);
    case 'arbitrary magnitude'
        %#function fdesign.arbmag
        constructor = 'arbmag';
    case 'arbitrary magnitude and phase'
        %#function fdesign.arbmagnphase
        constructor = 'arbmagnphase';
    case 'raised cosine'
        %#function fdesign.rcosine
        constructor = 'rcosine';
    case 'square root raised cosine'
        %#function fdesign.sqrtrcosine
        constructor = 'sqrtrcosine';
    otherwise
        % We still need to add the pragma for all new cases
        %#function fdesign.highpass
        %#function fdesign.lowpass
        %#function fdesign.bandpass
        %#function fdesign.bandstop
        %#function fdesign.nyquist
        %#function fdesign.differentiator
        %#function fdesign.hilbert
        %#function fdesign.halfband
        constructor = lower(dtype);
end

constructor = ['fdesign.' constructor];

% [EOF]
