classdef Resample < nirs.modules.AbstractModule
    %% RESAMPLE - Resamples time-series NIRS data.
    %
    % Options:
    %     Fs - new sampling frequency (Hz)
    
    
    properties
        Fs = 4; % new sampling frequency (Hz)
        Resample_Auxillary_Data = true;
    end
    properties (Hidden=true)
        antialias = [];
    end
    
    methods
        function obj = Resample( prevJob )
            obj.name = 'Resample';
            
            if nargin > 0
                obj.prevJob = prevJob;
            end
        end
        
        function data = runThis( obj, data )
            
            % Autodetect whether to use antialiasing (keeping consistent across
            % all files)
            if isempty(obj.antialias)
                obj.antialias = true;
                for i = 1:numel(data)
                    if length(data(i).time)>10000
                        obj.antialias = false;
                        break;
                    end
                end
            end
            
            for i = 1:numel(data)
                if(obj.Resample_Auxillary_Data && isa(data(i),'nirs.core.Data') && data(i).auxillary.count>0)
                    for j=1:data(i).auxillary.count
                        key= data(i).auxillary.keys{j};
                        try
                            data(i).auxillary(key)=obj.runThis( data(i).auxillary(key));
                        end
                    end
                end
                
                
                if obj.Fs < data(i).Fs
                    
                    % resample data
                    d = data(i).data;
                    t = data(i).time;
                    
                    % de-mean the data to avoid edge effects
                    mu = nanmean(d);
                    d = bsxfun(@minus,d,mu);
                    
                    if obj.antialias | verLessThan('signal','7.0')

                        % anti-aliasing filter
                        ord = floor( length(t) / 10 );
                        Fc = obj.Fs/data(i).Fs;
                        
                        b = fir1(ord,Fc);
                        
                        % zero phase filtering
                        d = filtfilt(b,1,d);
                        
                        % interpolation
                        N = floor((t(end)-t(1)) * obj.Fs)+1;
                        new_t = t(1) + (0:N-1)' / obj.Fs;
                        d = interp1(t,d,new_t,'linear','extrap');
                    else
                        % The data is too large, use the resample function
                        % instead
                        [d,new_t]=resample(d,t,obj.Fs);
                    end
                    
                    % restore original mean
                    d = bsxfun(@plus,d,mu);
                    
                    data(i).data = d;
                    data(i).time = new_t;
                end
            end
        end
    end
    
end

