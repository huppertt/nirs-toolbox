classdef Grangers < nirs.modules.AbstractModule
%% CONNECTIVITY - Computes all-to-all Multi-variate  Grangers connectivity model.
% Outputs nirs.core.ConnStats object

    properties
        modelorder;
    end
    methods
        function obj = Connectivity( prevJob )
           obj.name = 'Connectivity';
           obj.modelorder =[];
           
           if nargin > 0
               obj.prevJob = prevJob;
           end
        end
        
        function connStats = runThis( obj, data )
            
            for i = 1:numel(data)
                if(isempty(obj.modelorder))
                    pMax=data(i).Fs*4;
                else
                    pMax = obj.modelorder;
                end
                [G, F,dfe1,dfe2, P] = nirs.math.robust_mvgc(data(i).data,pMax);
                
                connStats(i)=nirs.core.ConnectivityStats();
                connStats(i).type = 'Correlation';
                connStats(i).description= ['Connectivity model of ' data(i).description];
                connStats(i).probe=data(i).probe;
                connStats(i).demographics=data(i).demographics;
                connStats(i).Grangers_F=F;
                connStats(i).dfe1=dfe1;
                connStats(i).dfe2=dfe2;
                connStats(i).Pearsons=[];
                
                disp(['Finished ' num2str(i) ' of ' num2str(length(data))]);
            end
        end
        
    end
    
end