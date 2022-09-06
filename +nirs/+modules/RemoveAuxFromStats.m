classdef RemoveAuxFromStats < nirs.modules.AbstractModule
    %% Remove Short Seperation - removes channels from probe/data
   
  
    
    methods
        function obj = RemoveAuxFromStats( prevJob )
            obj.name = 'RemoveAuxFromStats';
            
            if nargin > 0
                obj.prevJob = prevJob;
            end
        end
        
        function data = runThis( obj, data )
            for i = 1:numel(data)
                if(isa(data,'nirs.core.ChannelStats'))
                        lst=~cellfun('isempty', strfind(data(i).variables.cond,'Aux'));
                        data(i).variables(lst,:)=[];
                        data(i).beta(lst)=[];
                        data(i).covb(lst,:)=[];
                        data(i).covb(:,lst)=[];
                        
                    end
                end
            end
        end
    end
    

