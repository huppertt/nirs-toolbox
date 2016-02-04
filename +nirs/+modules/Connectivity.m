classdef Connectivity < nirs.modules.AbstractModule
%% CONNECTIVITY - Computes all-to-all connectivity model.
% Outputs nirs.core.ConnStats object

    properties
        corrfcn;  % function to use to compute correlation (see +nirs/+sFC for options)
    end
    methods
        function obj = Connectivity( prevJob )
           obj.name = 'Connectivity';
           obj.corrfcn = @(data)nirs.sFC.ar_corr(data,'4xFs',true);  %default to use AR-whitened robust correlation
           if nargin > 0
               obj.prevJob = prevJob;
           end
        end
        
        function connStats = runThis( obj, data )
            for i = 1:length(data)
                [r,p,dfe]=obj.corrfcn(data(i));
               
                connStats(i)=nirs.core.sFCStats();
                connStats(i).type = obj.corrfcn;
                connStats(i).description= ['Connectivity model of ' data(i).description];
                connStats(i).probe=data(i).probe;
                connStats(i).demographics=data(i).demographics;
     
                connStats(i).dfe=dfe;
                connStats(i).R=r;
                
                disp(['Finished ' num2str(i) ' of ' num2str(length(data))]);
            end
        end
        
    end
    
end
