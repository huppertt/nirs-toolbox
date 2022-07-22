classdef ApplyStatsContrast < nirs.modules.AbstractModule
    %% ApplyStatsContrast
    %   This function applies a t-test contrast as part of a job pipeline
    %
    % Options:
    %     Contrasts- a cell array of contrast vectors, contrast strings, or
    %     function handles.  Function handles must return a valid contrast
    %     vector (e.g. nirs.design.contrastvector() using StatsData as input)
    
    
    properties
        Contrasts={};
    end
    
    methods
        function obj = ApplyStatsContrast( prevJob )
            obj.name = 'ApplyStatsContrast';
            
            if nargin > 0
                obj.prevJob = prevJob;
            end
        end
        
        function data = runThis( obj, data )
            for id=1:length(data)
                
                if(~isempty(obj.Contrasts))
%                     C=zeros(0,length(data(id).conditions));
%                     for i=1:length(obj.Contrasts)
%                         if(isstr(obj.Contrasts{i}) | iscellstr(obj.Contrasts{i}))
%                             C=[C; nirs.design.contrastvector(obj.Contrasts{i},data(id).conditions)];
%                         elseif(isa(obj.Contrasts{i},'function_handle'))
%                             a=obj.Contrasts{i}(data(id));
%                             if(isstr(a) | iscellstr(a))
%                                 C=[C; nirs.design.contrastvector(a,data(id).conditions)];
%                             else
%                                 C=[C; a];
%                             end
%                         else
%                             C=[C; obj.Contrasts{i}];
%                         end
%                     end
                    
                    
                    data(id) = data(id).ttest(obj.Contrasts);
                end
            end
        end
        
    end
    
end

