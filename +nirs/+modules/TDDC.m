classdef TDDC < nirs.modules.AbstractModule
%% Temporal Derivative Distribution Clamping - Attempts to shrink artifactual shifts and spikes.
% 
% Options: 
%     lowpass_cutoff - Filter cutoff in Hz for separating high and low
%     frequencies (only low is corrected, then they are recombined)

    properties
        lowpass_cutoff = 0.5; % Filter cutoff to prevent high-frequency noise from inflating the stddev
    end
    
    methods
        function obj = TDDC( prevJob )
           obj.name = 'Temporal Derivative Distribution Clamping';
           
           if nargin > 0
               obj.prevJob = prevJob;
           end
        end
        
        function data = runThis( obj, data )
            
            for i = 1:length(data)
            
                % Separate into high and low frequency components
                Wn = obj.lowpass_cutoff * 2/data(i).Fs;
                [fb,fa]=butter( 3 , Wn );
                low = filtfilt( fb , fa , data(i).data );
                high = data(i).data - low;
                
                % Perform TDDC on low component
                low = nirs.math.tddc( low );
                
                % DCT detrending of corrected low component
                k  = size(data(i).data, 1); 
                RT = 1/data(i).Fs; 
                n  = fix(2*(k*RT)/128 + 1);
                X0 = spm_dctmtx(k,n);
                X0 = X0(:,2:end);
                low = low - X0 * (X0' * low);
                
                % Merge corrected low w/ original high, restoring original
                % mean
                low = bsxfun( @minus , low , mean(low) );
                high = bsxfun( @minus , high , mean(high) );
                data(i).data = bsxfun( @plus , low+high , mean(data(i).data) );
                
            end
            
        end
        
    end
    
end