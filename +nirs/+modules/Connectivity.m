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
                pMax=round(4*data(i).Fs);
                
              [G, F, df1, df2, p] = nirs.math.mvgc(data(i).data,pMax);
                %G=log(sqrt(F.*df1./df2+1))

                
                [yf,a]= nirs.math.innovations(data(i).data,pMax);
                [rcorr,pcorr]=corrcoef(yf);
                
                connStats(i)=nirs.core.ConnectivityStats();
                connStats(i).description= ['Connectivity model of ' data(i).description];
                connStats(i).probe=data(i).probe;
                connStats(i).demographics=data(i).demographics;
                connStats(i).Grangers_F=F;
                connStats(i).dfe1=df1;
                connStats(i).dfe2=df2;
                connStats(i).Pearsons=rcorr;
                
                disp(['Finished ' num2str(i) ' of ' num2str(length(data))]);
            end
        end
        
    end
    
end
