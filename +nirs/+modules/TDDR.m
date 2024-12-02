classdef TDDR < nirs.modules.AbstractModule
%% Temporal Derivative Distribution Repair - Corrects motion artifacts by downweighting outlier fluctuations

properties    
    usePCA = false;  % Do correction on the PrinComp of the data instead of channel space
    split_PosNeg=false;  % Process pos and negative derivatives seperately
end
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

                    if(obj.usePCA)
                        [U,S,V]=nirs.math.mysvd(data(i).data);
                        U=nirs.math.tddr( U , data(i).Fs,obj.split_PosNeg );
                        data(i).data=U*S*V';
                        
                    else
                        data(i).data = nirs.math.tddr( data(i).data , data(i).Fs ,obj.split_PosNeg);
                    end
                catch
                    disp(['Error in file ' num2str(i) ' ' lasterr]);
                end
            end
            
        end
        
    end
    
end
