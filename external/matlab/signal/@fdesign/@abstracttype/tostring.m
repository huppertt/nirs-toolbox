function str = tostring(this)
%TOSTRING

%   Copyright 2005-2013 The MathWorks, Inc.

s = get(this);
s = reorderstructure(s, 'Response', 'Specification', 'Description');

% Remove Description field
s = rmfield(s, 'Description');

% Remove sampling frequencies and normalized
s = rmfield(s, 'Fs');
if isfield(s,'Fs_in'), s = rmfield(s, 'Fs_in'); end
if isfield(s,'Fs_out'), s = rmfield(s, 'Fs_out'); end

% Remove CICRateChangeFactor if it exists and its value is 1
if isfield(s,'CICRateChangeFactor') && this.CICRateChangeFactor == 1
  s = rmfield(s,'CICRateChangeFactor');
end

s = rmfield(s,'NormalizedFrequency');

specs = getcurrentspecs(this);
if ~isempty(specs) && isfromdesignfilt(specs)
  s = rmfield(s,'Specification');
end

f = fieldnames(s);

param = cell(length(f)+1, 1);
value = param;

param{1} = getString(message('signal:dfilt:info:SamplingFrequency'));
if this.NormalizedFrequency
    if ~isempty(specs) && isfromdesignfilt(specs)
      value{1} = ['2 (' getString(message('signal:dfilt:info:Normalized')) ')'] ;
    else
      value{1} = getString(message('signal:dfilt:info:NAnormalizedFrequency'));
    end
    m        = 1;
else
    
    [fs, m, prefix] = engunits(this.Fs);
    
    fpost    = sprintf(' %sHz', prefix);
    value{1} = [num2str(fs) fpost];
end

for indx = 1:length(f)
    
    if isnumeric(s.(f{indx})),
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
            case {'fcutoff', 'fcutoff1', 'fcutoff2'}
                d = [d getString(message('signal:dfilt:info:cutoffFrequency'))];                
            case {'groupdelay'}
                d = [d getString(message('signal:dfilt:info:GroupDelay'))];
            case {'nomgrpdelay'}
                d = [d getString(message('signal:dfilt:info:NominalGroupDelay'))];                
            case {'filterorder'}
                d = [d getString(message('signal:dfilt:info:FilterOrder'))];       
            case {'nbands'}
                d = [d getString(message('signal:dfilt:info:NumberofBands'))];
            case {'cicratechangefactor'}            
                  d = [d getString(message('signal:dfilt:info:CICRateChangeFactor'))];                
            case {'numberofsections'}            
                  d = [d getString(message('signal:dfilt:info:NumberofSections'))];    
            case {'differentialdelay'}
                  d = [d getString(message('signal:dfilt:info:DifferentialDelay'))];
            case {'decimationfactor'}
                  d = [d getString(message('signal:dfilt:info:DecimationFactor'))];    
            case {'interpolationfactor'}
                  d = [d getString(message('signal:dfilt:info:InterpolationFactor'))];                      
            otherwise
                  d = getTranslatedString('signal:dfilt:info',f{indx});
        end
        
        BiAmplitudesPattern = 'B.Amplitudes';
        BiRipplePattern = 'B.Ripple';
        BiFreqRespPattern = 'B.FreqResponse';
        if any(strcmpi(f{indx},{'Amplitudes','FreqResponse','FilterOrder',...
                'NumOrder','DenOrder','NumberOfSections','DifferentialDelay',...
                'DecimationFactor','InterpolationFactor','q','BandsPerOctave',...
                'NBands','Ripple','Class','CICRateChangeFactor'})) || any(regexp(f{indx},BiAmplitudesPattern)) ||...
                any(regexp(f{indx},BiRipplePattern)) || any(regexp(f{indx},BiFreqRespPattern))    
            post = '';
            value{indx+1} = num2str(this.(f{indx}));
        elseif strncmpi(f{indx}, 'a', 1) || (strncmpi(f{indx}, 'g', 1) && ~strcmpi(f{indx},'groupdelay'))
            post = ' dB';
            value{indx+1} = num2str(this.(f{indx}));
        elseif any(strfind(lower(f{indx}), 'delay'))
            delay = this.(f{indx});
            post    = sprintf(' %s', getString(message('signal:dfilt:info:samples')));
            value{indx+1} = num2str(delay);
            
            if ~this.NormalizedFrequency,
              [~, m, prefix] = engunits(delay,'time');
              post = sprintf(' %s', prefix);
            end
            value{indx+1} = num2str(delay*m);
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
        
        param{indx+1} = d;
    elseif islogical(s.(f{indx}))
        param{indx+1} = f{indx};
        if s.(f{indx})
          value{indx+1} = getString(message('signal:dfilt:info:true'));
        else
          value{indx+1} = getString(message('signal:dfilt:info:false'));
        end          
    else
        if strcmpi(f{indx},'weightingtype')
          param{indx+1} = getString(message('signal:dfilt:info:WeightingType'));
        elseif strcmpi(f{indx},'multiratetype')
          param{indx+1} = getString(message('signal:dfilt:info:MultirateType'));          
        else
          param{indx+1} = getTranslatedString('signal:dfilt:info',f{indx});
        end
        value{indx+1} = s.(f{indx});
    end
end

str = [char(param) repmat(' : ', length(param), 1) char(value)];


% [EOF]

