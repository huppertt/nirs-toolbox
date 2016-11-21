classdef TrimBaseline < nirs.modules.AbstractModule
%% TrimBaseline - Removes excessive baseline at the beginning or end of file.
%
% Options:
%     preBaseline  - maximum baseline (seconds) at the beginning of scan
%     postBaseline - maximum baseline (seconds) after final task period ends
    
    properties
        preBaseline  = 30;  % maximum baseline (seconds) at the beginning of scan
        postBaseline = 30;  % maximum baseline (seconds) after final task period ends
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
                
%                 s = sum( abs(s), 2 );
%                 
%                 % find first/last stim period and calculate time inverval
%                 t_min = t( find(s>0,1,'first') ) - obj.preBaseline;
%                 t_max = t( find(s>0,1,'last' ) ) + obj.postBaseline;
%                 

                t_min=min(vertcat(stims{:}.onset))- obj.preBaseline;
                t_max=max(vertcat(stims{:}.onset)+vertcat(stims{:}.dur))+ obj.postBaseline;

                % extract data from the time inverval t_min to t_max
                if(~isempty(t_min))
                lst = t >= t_min & t <= t_max;
                t = t(lst);
                d = d(lst,:);
                end
                
                data(i).data = d;
                data(i).time = t;
                
            end
        end
    end
    
end

