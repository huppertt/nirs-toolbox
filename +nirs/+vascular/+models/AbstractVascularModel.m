classdef AbstractVascularModel
    % Abstract class for vascular models
    
    properties
        parameters
        outputNames;
        description;
        
        fitter;
    end
     properties(Access=protected)
        states;
        name;
    end
    methods 
        function obj = AbstractVascularModel
            obj.outputNames = {'CMRO2','CBF'};
            obj.fitter=@nirs.vascular.estimation.extended_kalman_fit;
        end
        function varargout = StateEstimate(obj,data)
            % This function estimates the model
            % states
            % This returns the states in obj.outputNames
            states=obj.fitter(obj.model,data);
            
            for i=1:length(obj.outputNames)
                id=find(ismember({states.name},obj.outputNames{i}));
            
                varargout{i}=states(id).data;
            end
            
            
        end
        
        function data = sim(obj,input)
            nlgr=obj.model;
            opt = simOptions('InitialCondition',cell2mat(getinit(nlgr)));
            data=sim(nlgr,input,opt); 
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