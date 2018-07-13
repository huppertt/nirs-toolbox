classdef RemoveShortSeperationRegressors < nirs.modules.AbstractModule
    %% AddShortSeperationRegressors - Adds short seperation data as regressors to the GLM model
    %
    
   methods
        function obj = RemoveShortSeperationRegressors( prevJob )
            obj.name = 'RemoveShortSeperationRegressors';
           if nargin > 0
                obj.prevJob = prevJob;
            end
        end
        
        function data = runThis( obj, data )
            
            if(isa(data,'nirs.core.ChannelStats'))
            for i = 1:numel(data)
                names=nirs.getStimNames(data(i));
                names2={}; cnt=1;
                for j=1:length(names)
                    if(isempty(strfind(names{j},'SS_PCA')))
                        names2{cnt}=names{j};
                        cnt=cnt+1;
                    end
                end
               data(i)=data(i).ttest(names2);
            end
            % this removes the SS data/ttests from the stats 
            j=nirs.modules.RemoveShortSeperations;
            data=j.run(data);
            else
                for i = 1:numel(data)
                    names=nirs.getStimNames(data(i));
                    names2={}; cnt=1;
                    for j=1:length(names)
                        if(isempty(strfind(names{j},'SS_PCA')))
                            names2{cnt}=names{j};
                            cnt=cnt+1;
                        end
                    end
                    j=nirs.modules.KeepStims;
                    j.listOfStims=names2;
                    data(i)=j.run(data(i));
                end
            end
        end
    end
    
end

