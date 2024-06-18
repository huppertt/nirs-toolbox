classdef Resize_Heads < nirs.modules.AbstractModule
%% Resize the probe for each subject based on the head size stored in the demogrpahics
%
% Usage:
% job = nirs.modules.Resize_Heads();

    properties
        demographics_key='circumference';
    end
    
    methods
        function obj = Resize_Heads( prevJob )
         	if nargin > 0, obj.prevJob = prevJob; end
            obj.name = 'Resize Heads';
        end
        
        function data = runThis( obj , data )
           
            for i = 1:length(data)
                headsize_tmp = data(i).demographics(obj.demographics_key);
                if(isa(headsize_tmp,'Dictionary'))
                    headsize=headsize_tmp;
                elseif(isa(headsize_tmp,'struct'))
                    headsize=Dictionary(headsize_tmp);
                else
                    headsize=Dictionary({'circumference'},{headsize_tmp});
                end


                data(i).probe = nirs.util.resize_fixedprobe(data(i).probe,headsize);
            end
            
        end
        
    end
end