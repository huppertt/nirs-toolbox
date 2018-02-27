classdef Resample < nirs.modules.AbstractModule
    %% RESAMPLE - Resamples time-series NIRS data.
    %
    % Options:
    %     Fs - new sampling frequency (Hz)
    
    
    properties
        Fs = 4; % new sampling frequency (Hz)
    end
    
    methods
        function obj = Resample( prevJob )
            obj.name = 'Resample';
            
            if nargin > 0
                obj.prevJob = prevJob;
            end
        end
        
        function data = runThis( obj, data )
            for i = 1:numel(data)
                if obj.Fs < data(i).Fs
                    
                    % resample data
                    d = data(i).data;
                    t = data(i).time;
                    N = floor((t(end)-t(1)) * obj.Fs)+1;
                    new_t = t(1) + (0:N-1)' / obj.Fs;
                    
                    % de-mean the data to avoid edge effects
                    mu = nanmean(d);
                    d = bsxfun(@minus,d,mu);
                    
                    if(length(t)<10000)

                        % anti-aliasing filter
                        ord = floor( length(t) / 10 );
                        Fc = obj.Fs/data(i).Fs;
                        
                        b = fir1(ord,Fc);
                        
                        % zero phase filtering
                        d = filtfilt(b,1,d);
                        
                        % interpolation
                        d = interp1(t,d,new_t,'linear');
                    else
                        % The data is too large, use the resample function
                        % instead
                        [P,Q] = rat(obj.Fs/data(i).Fs);
                        d=resample(d,P,Q);
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

