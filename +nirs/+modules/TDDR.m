classdef TDDR < nirs.modules.AbstractModule
%% Temporal Derivative Distribution Repair - Corrects motion artifacts by downweighting outlier fluctuations
    methods
        function obj = TDDR( prevJob )
           obj.name = 'Temporal Derivative Distribution Repair';
           obj.citation=['Fishburn, Frank A., Ludlum, Ruth S., Vaidya, Chandan J., and Medvedev, Andrei V. '...
                '"Temporal Derivative Distribution Repair (TDDR): A motion correction method for fNIRS." '...
                'NeuroImage 184 (2019): 171-179.'];
           if nargin > 0
               obj.prevJob = prevJob;
           end
        end
        
        function data = runThis( obj, data )
            
            for i = 1:numel(data)
                try
                    data(i).data = nirs.math.tddr( data(i).data , data(i).Fs );
                catch
                    disp(['Error in file ' num2str(i) ' ' lasterr]);
                end
            end
            
        end
        
    end
    
end
