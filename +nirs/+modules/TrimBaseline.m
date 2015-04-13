classdef TrimBaseline < nirs.modules.AbstractModule
% Trim pre/post baseline
    
    properties
        preBaseline  = 30;
        postBaseline = 30;
    end
    
    methods
        function obj = TrimBaseline( prevJob )
           obj.name = 'Trim Pre/Post Baseline';
           if nargin > 0
               obj.prevJob = prevJob;
           end
        end
        
        function data = runThis( obj, data )
            for i = 1:length(data)
                
                d = data(i).data;
                t = data(i).time;
                
                % get stims 
                stims = data(i).stimulus.values;
                
                % get stim vectors from stims
                s = zeros(size(t,1),length(stims));
                for j = 1:length(stims)
                    s(:,j) = stims{j}.getStimVector( t );
                end
                
                s = abs( sum(s,2) );
                
                % find first/last stim period and calculate time inverval
                t_min = t( find(s>0,1,'first') ) - obj.preBaseline;
                t_max = t( find(s>0,1,'last' ) ) + obj.postBaseline;
                
                % extract data from the time inverval t_min to t_max
                lst = t >= t_min & t <= t_max;
                t = t(lst);
                d = d(lst,:);
                
                data(i).data = d;
                data(i).time = t;
                
            end
        end
    end
    
end

