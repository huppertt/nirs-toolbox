classdef Connectivity < nirs.modules.AbstractModule
%% CONNECTIVITY - Computes all-to-all connectivity model.
% Outputs nirs.core.ConnStats object

    
    methods
        function obj = Connectivity( prevJob )
           obj.name = 'Connectivity';
          
           if nargin > 0
               obj.prevJob = prevJob;
           end
        end
        
        function connStats = runThis( obj, data )
            for i = 1:length(data)
                pMax=data(i).Fs*4;
                [yf,a]= nirs.math.innovations(data(i).data,pMax);
                [rcorr,pcorr]=nirs.math.robust_corrcoef(yf);
                dfe2=repmat(size(yf,1)-1,1,size(yf,2));
                
                connStats(i)=nirs.core.ConnectivityStats();
                connStats(i).type = 'Correlation';
                connStats(i).description= ['Connectivity model of ' data(i).description];
                connStats(i).probe=data(i).probe;
                connStats(i).demographics=data(i).demographics;
                connStats(i).Grangers_F=[];
                connStats(i).dfe1=[];
                connStats(i).dfe2=dfe2;
                connStats(i).Pearsons=rcorr;
                
                disp(['Finished ' num2str(i) ' of ' num2str(length(data))]);
            end
        end
        
    end
    
end
