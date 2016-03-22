classdef AbstractVascularModel
    % Abstract class for vascular models
    
    properties
        parameters
        outputNames;
        description;
    end
     properties(Access=protected)
        model;
        states;
        fitfunction;
        name;
    end
    methods 
        function obj = AbstractVascularModel
            obj.outputNames = {'CMRO2','CBF'};
        end
        function varargout = StateEstimate(obj,data)
            % This function estimates the model
            % states
            % This returns the states in obj.outputNames
                        
        end
        function obj=set.outputNames(obj,Names)
            % Change the output of the model and make sure it is valid
          if(isempty(obj.states))
                obj.outputNames=Names;
           else
               lst=find(~ismember(Names,obj.states));
               if(~isempty(lst))
                   warning('Requested name is not a state');
                   return
               end
               obj.outputNames=Names;
           end
        end
    end
end