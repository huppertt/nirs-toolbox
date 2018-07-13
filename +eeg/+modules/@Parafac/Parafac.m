classdef Parafac < nirs.modules.AbstractModule

    properties
        maxcomponents;
        datacombine;
    end
    
    methods

        function obj =Parafac( prevJob )
           obj.name = 'Parallel factor analysis model';
           if nargin > 0
               obj.prevJob = prevJob;
           end
           obj.maxcomponents=4;
           obj.datacombine=false;
        end
        
        function data = runThis( obj, data )
            if(~obj.datacombine)
                for idx=1:length(data)
                    data(idx)=obj.computefactors(data(idx));
                end
                return;
            else
                data=obj.computefactors(data);
            end
            
        end
        data=computefactors(obj,data);
    end
end