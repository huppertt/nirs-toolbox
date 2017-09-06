classdef variablecontainer
    properties
        type
        value
    end
    methods
        function obj=variablecontainer(obj,val)
            obj.type=class(val);
            obj.value=val;
        end
        function out=get(obj)
            out=obj.value;
        end
        function set(obj,val)
            obj.value=val;
        end
    end
end