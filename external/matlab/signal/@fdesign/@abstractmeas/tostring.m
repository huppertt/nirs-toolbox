function str = tostring(this)
%TOSTRING

%   Author(s): J. Schickler
%   Copyright 2005-2010 The MathWorks, Inc.

% Get all of the field names
f = fieldnames(this);

% Remove 'NormalizedFrequency and Fs because these are special cased above.
f(1:2) = [];

param = cell(length(f)+1, 1);
value = param;

param{1} = getString(message('signal:dfilt:info:SamplingFrequency'));
if this.NormalizedFrequency
    value{1} = getString(message('signal:dfilt:info:NAnormalizedFrequency'));
    fpost    = '';
    m        = 1;
else
    
    [fs m prefix] = engunits(this.Fs);
    
    fpost    = sprintf(' %sHz', prefix);
    value{1} = [num2str(fs) fpost];
end

for indx = 1:length(f)

    num = str2num(f{indx}(end)); %#ok<ST2NM>
    if isempty(num),
        d = '';
    else
        switch num
            case 1
                d = [getString(message('signal:dfilt:info:First')) ' '];
            case 2
                d = [getString(message('signal:dfilt:info:Second')) ' '];
            case 3
                d = [getString(message('signal:dfilt:info:Third')) ' '];
            case 4
                d = [getString(message('signal:dfilt:info:Fourth')) ' '];
        end
    end

    switch lower(f{indx})
        case {'transitionwidth', 'transitionwidth1', 'transitionwidth2'};
            d    = [d getString(message('signal:dfilt:info:TransitionWidth'))]; %#ok<*AGROW>
            post = fpost; %#ok<*NASGU>
        case {'lowtransitionwidth'};
            d    = [d getString(message('signal:dfilt:info:LowTransitionWidth'))];
            post = fpost;  
        case {'hightransitionwidth'};
            d    = [d getString(message('signal:dfilt:info:HighTransitionWidth'))];
            post = fpost;                        
        case {'apass', 'apass1', 'apass2'}
            d = [d getString(message('signal:dfilt:info:PassbandRipple'))];
        case {'astop', 'astop1', 'astop2'}
            d = [d getString(message('signal:dfilt:info:StopbandAtten'))];
        case {'fpass', 'fpass1', 'fpass2'}
            d = [d getString(message('signal:dfilt:info:PassbandEdge'))];
        case {'fstop', 'fstop1', 'fstop2'}
            d = [d getString(message('signal:dfilt:info:StopbandEdge'))];
        case {'f3db', 'f6db', 'f3db1', 'f6db1', 'f3db2', 'f6db2'}
            d = [d f{indx}(2) getString(message('signal:dfilt:info:dBPoint'))];
        case {'nomgrpdelay'}
            d = getString(message('signal:dfilt:info:NominalGroupDelay'));
        case {'fracdelayerror'}
            d = getString(message('signal:dfilt:info:FractionalDelayError'));
        case {'f0'},
            d = getString(message('signal:dfilt:info:CenterFrequency'));
        case {'bw'},
            d = getString(message('signal:dfilt:info:Bandwidth'));
        case {'g0'},
            d = getString(message('signal:dfilt:info:CenterFrequencyGain'));
        case {'gref'},
            d = getString(message('signal:dfilt:info:ReferenceGain'));
        case {'gbw'},
            d = getString(message('signal:dfilt:info:BandwidthGain'));
        case {'gpass'},
            d = getString(message('signal:dfilt:info:PassbandGain'));
        case {'gstop'},
            d = getString(message('signal:dfilt:info:StopbandGain'));
        case {'bwpass'},
            d = getString(message('signal:dfilt:info:PassbandBandwidth'));
        case {'bwstop'},
            d = getString(message('signal:dfilt:info:StopbandBandwidth'));
        case {'notchfrequencies'},
            d = getString(message('signal:dfilt:info:NotchFrequencies'));
         case {'peakfrequencies'},
            d = getString(message('signal:dfilt:info:PeakFrequencies'));
        case {'gnotch'},
            d = getString(message('signal:dfilt:info:NotchGains'));
         case {'gpeak'},
            d = getString(message('signal:dfilt:info:PeakGains'));
         case {'qa'},
            d = getString(message('signal:dfilt:info:QualityFactoraudio'));            
        case {'totalgroupdelay'}
            d = getString(message('signal:dfilt:info:TotalGroupDelay')); 
        case {'freqresponse'}
            d = getString(message('signal:dfilt:info:FrequencyResponse'));  
        otherwise
            d = getTranslatedString('signal:dfilt:info',f{indx});
    end
    if isempty(this.(f{indx}))
        if strcmpi(f{indx},'qa')
            value{indx+1} = getString(message('signal:dfilt:info:NotMeasurableForN2'));            
        else
            value{indx+1} = getString(message('signal:dfilt:info:Unknown'));
        end
    else
        if strcmpi(f{indx}, 'magnitudes') || strncmpi(f{indx}, 'g', 1) || (strncmpi(f{indx}, 'a', 1) && ~strcmpi(f{indx}, 'amplitudes'))
            post = ' dB';
            value{indx+1} = num2str(this.(f{indx}));
        elseif any(strfind(lower(f{indx}), 'delay'))
            delay = this.(f{indx});
            post    = sprintf(' %s', getString(message('signal:dfilt:info:samples')));            
            value{indx+1} = num2str(delay);  
            
            if exist('fs') %#ok<EXIST>
              [~, m prefix] = engunits(delay,'time');
              post = sprintf(' %s', prefix);
            end
            value{indx+1} = num2str(delay*m);
        elseif any(strfind(lower(f{indx}), 'amplitudes')) || any(strfind(lower(f{indx}), 'freqresponse')) || strcmpi(f{indx},'q') || strcmpi(f{indx},'qa')
            value{indx+1} = num2str(this.(f{indx}));
            post    = '';
        else
            if this.NormalizedFrequency
                value{indx+1} = num2str(this.(f{indx}));
                post    = '';
            else
                [val , ~, valprefix] = engunits(this.(f{indx}));
                post    = sprintf(' %sHz', valprefix);
                value{indx+1} = num2str(val);
            end
        end

        value{indx+1} = sprintf('%s%s', value{indx+1}, post);
    end
    param{indx+1} = d;
end

str = [char(param) repmat(' : ', length(param), 1) char(value)];

% [EOF]
