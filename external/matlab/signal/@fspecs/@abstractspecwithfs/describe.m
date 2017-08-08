function p = describe(this)
%DESCRIBE   Describe the object.

%   Author(s): J. Schickler
%   Copyright 1988-2008 The MathWorks, Inc.

p      = propstoadd(this);
p(1:2) = [];

for indx = 1:length(p),
    
    % Special case FilterOrder and TransitionWidth.
    switch lower(p{indx})
        case 'bwpass'
            d = 'Passband Width';
        case 'bwstop'
            d = 'Stopband Width';
        case 'polyphaselength'
            d = 'Polyphase Length';
        case 'filterorder',
            d = 'Filter Order';
        case 'polynomialorder',
            d = 'Polynomial Order';            
        case 'transitionwidth';
            d = 'Transition Width';
        case 'f0'
            d = 'Center Frequency';
        case 'fc'
            d = 'Cutoff Frequency';            
        case 'bw'
            d = 'Bandwidth at GBW';
        case 'gref'
            d = 'Reference Gain (dB)';
        case 'g0'
            d = 'Gain at Center Frequency (dB)';
        case 'gbw'
            d = 'Bandwidth Gain (dB)';
        case 'gpass'
            d = 'Passband Gain (dB)';
        case 'gstop'
            d = 'Stopband Gain (dB)';
        case 'q'
            d = 'Quality Factor';
        case 'qa'
            d = 'Quality Factor (audio)';
        case 's'
            d = 'Shelf Slope Parameter';            
        case 'rollofffactor'
            d = 'Rolloff Factor';
        case 'numberofsymbols'
            d = 'Filter Order in Symbols';
        case 'bt'
            d = 'Bandwidth-time product';
        otherwise
            num = str2double(p{indx}(end));
            if isnan(num),
                d = '';
            else
                p{indx}(end) = [];
                switch num
                    case 1
                        d = 'First ';
                    case 2
                        d = 'Second ';
                    case 3
                        d = 'Third ';
                    case 4
                        d = 'Fourth ';
                end
            end
            switch p{indx}(2:end)
                case 'pass'
                    ispass = true;
                    d = sprintf('%sPassband ', d);
                case 'stop'
                    ispass = false;
                    d = sprintf('%sStopband ', d);
                case 'cutoff'
                    ispass = false;
                    d = sprintf('%sCutoff ', d);
                case '3dB'
                    ispass = false;
                    d = sprintf('%s3dB ', d);
            end
            switch p{indx}(1)
                case 'A'
                    if ispass
                        d = sprintf('%sRipple (dB)', d);
                    else
                        d = sprintf('%sAttenuation (dB)', d);
                    end
                case 'F'
                    d = sprintf('%sFrequency', d);
                case 'W'
                    d = sprintf('%sWeight', d);
            end
    end
    p{indx} = d;
end

% Make sure that we get a column.
p = p(:);

% [EOF]
