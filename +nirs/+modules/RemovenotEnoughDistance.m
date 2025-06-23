classdef RemovenotEnoughDistance < nirs.modules.AbstractModule
    %% Remove not Enough Distance - removes channels from probe/data
   
  
    
    methods
        function obj = RemovenotEnoughDistance( prevJob )
            obj.name = 'RemovenotEnoughDistance';
            
            if nargin > 0
                obj.prevJob = prevJob;
            end
        end
        
        function data = runThis( obj, data )
            for i = 1:numel(data)
                if(ismember('notEnoughDistance',data(i).probe.link.Properties.VariableNames));
                    if(isa(data,'nirs.core.Data'))
                        lst=find(data(i).probe.link.notEnoughDistance);
                        data(i).probe.link(lst,:)=[];
                        data(i).data(:,lst)=[];
                    elseif(isa(data,'nirs.core.ChannelStats'))
                        lst=find(data(i).probe.link.notEnoughDistance);
                        data(i).probe.link(lst,:)=[];
                        lst=find(data(i).variables.notEnoughDistance);
                        data(i).variables(lst,:)=[];
                        data(i).beta(lst)=[];
                        data(i).covb(lst,:)=[];
                        data(i).covb(:,lst)=[];
                        
                    end
                end
            end
        end
    end
    
end

