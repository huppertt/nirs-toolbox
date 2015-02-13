classdef Resample < nirs.functional.AbstractModule
    %UNTITLED3 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        Fs = 4;
    end
    
    methods

        function obj = Resample( prevJob )
           obj.name = 'Resample';
           
           if nargin > 0
               obj.prevJob = prevJob;
           end
        end
        
        function data = execute( obj, data )
            for i = 1:length(data)
                
                assert( obj.Fs < data(i).Fs )
                
                % resample data
                d = data(i).data;
                t = data(i).time;
                
                N = floor(t(end) * obj.Fs);
                new_t = t(1) + (0:N-1)' / obj.Fs;
                
                % anti-aliasing filter
                ord = floor( length(new_t) / 10 );
                Fc = obj.Fs/data(i).Fs;
                [zz,pp,kk] = butter(ord,Fc,'low');  % Butterworth filter
                [sos,g] = zp2sos(zz,pp,kk);          % Convert to SOS form
                
                d = filtfilt(sos,g,d);
                
                % interpolation
                d = interp1(t,d,new_t);

                data(i).data = d;
                data(i).time = new_t;
                
%                 % resample stims
%                 stimulus = data(i).stimulus;
%                	keys = stimulus.keys;
%                 
%                 for j = 1:length(keys)
%                     %stimulus( keys(j) ) = filtfilt( sos,g,stimulus(keys(j)) );
%                     stimulus( keys{j} ) = interp1(t, stimulus(keys{j}), new_t); 
%                 end
%                 data(i).stimulus = stimulus;
                
            end
        end
        
        function options = getOptions( obj )
            options = [];
        end
           
        function obj = putOptions( obj, options )
        end
        
    end
    
end

