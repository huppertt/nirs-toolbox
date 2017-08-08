function nfft = calcnfft(hopts,segLen)
%CALCNFFT   Get the numeric value of NFFT even when it's set to a string.

%   Copyright 1988-2006 The MathWorks, Inc.

freqpoints = hopts.FreqPoints;
switch lower(freqpoints),
    case 'all',
        nfft = hopts.NFFT;
        if ischar(nfft),
            switch lower(nfft),
                case 'nextpow2',
                    nfft = max(256,2^nextpow2(segLen)); % segLen=input length for all but welch.
                case 'auto',
                    nfft = max(256,segLen);
            end
        end
    case 'user defined',
        nfft = hopts.frequencyvector;
        if ischar(nfft),
            switch lower(nfft),
                case 'auto',
                    nfft = max(256,segLen);
            end
            Fs = hopts.Fs;
            if ischar(Fs)
                if strcmp(Fs, 'Normalized')
                    Fs = 2*pi;
                end
            end

            range = 'whole';
            if ishalfnyqinterval(hopts),
                range = 'half';
            end
            nfft = psdfreqvec('npts',nfft, 'Fs',Fs, 'Range',range);
        end
end
% [EOF]
