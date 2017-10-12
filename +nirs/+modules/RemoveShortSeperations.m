classdef RemoveShortSeperation < nirs.modules.AbstractModule
    %% Remove Short Seperation - removes channels from probe/data
   
  
    
    methods
        function obj = RemovelShortSeperation( prevJob )
            obj.name = 'RemoveShortSeperation';
            
            if nargin > 0
                obj.prevJob = prevJob;
            end
        end
        
        function data = runThis( obj, data )
            for i = 1:numel(data)
              
               if(isa(data,'nirs.core.Data'))
                    lst=find(data(i).probe.link.ShortSeperation);
                    data(i).probe.link(lst,:)=[];
                    data(i).data(:,lst)=[];
               elseif(isa(data,'nirs.core.ChannelStats'))
                    lst=find(data(i).probe.link.ShortSeperation);
                    data(i).probe.link(lst,:)=[];
                    lst=find(data(i).variables.ShortSeperation);
                    data(i).variables(lst,:)=[];
                    data(i).beta(lst)=[];
                    data(i).covb(lst,:)=[];
                    data(i).covb(:,lst)=[];
                    
               end
            end
        end
    end
    
end

