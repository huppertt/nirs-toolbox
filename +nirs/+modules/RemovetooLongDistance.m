classdef RemovetooLongDistance < nirs.modules.AbstractModule
    %% Remove too long Distance - removes channels from probe/data
   
  
    
    methods
        function obj = RemovetooLongDistance( prevJob )
            obj.name = 'RemovetooLongDistance';
            
            if nargin > 0
                obj.prevJob = prevJob;
            end
        end
        
        function data = runThis( obj, data )
            for i = 1:numel(data)
                if(ismember('tooLongDistance',data(i).probe.link.Properties.VariableNames));
                    if(isa(data,'nirs.core.Data'))
                        lst=find(data(i).probe.link.tooLongDistance);
                        data(i).probe.link(lst,:)=[];
                        data(i).data(:,lst)=[];
                    elseif(isa(data,'nirs.core.ChannelStats'))
                        lst=find(data(i).probe.link.tooLongDistance);
                        data(i).probe.link(lst,:)=[];
                        lst=find(data(i).variables.tooLongDistance);
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

