classdef AR_Prewhiten < nirs.modules.AbstractModule
%% AR_Prewhiten - Remove serial correlations from data
%
% Options:
%     modelorder - % Maximum model order (in seconds)

    properties
         modelorder = 4; % Maximum model order
         verbose = false;
    end
    
    methods

        function obj =  AR_Prewhiten( prevJob )
           obj.name = 'Autoregressive prewhitening';         
           if nargin > 0
               obj.prevJob = prevJob;
           end
        end
        
        function data = runThis( obj, data )
            
            nfiles = length(data);
            
            parfor i = 1:nfiles
                
                if obj.verbose
                    fprintf('Prewhitening file %i of %i...\n',i,nfiles);
                end
                
                p = data(i).Fs * obj.modelorder;
                data(i).data = nirs.math.innovations( data(i).data , p );
                
            end
            
        end
    end
    
end

