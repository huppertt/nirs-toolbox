classdef TrimBaseline < nirs.functional.AbstractModule
    %UNTITLED6 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        preBaseline  = 40;
        postBaseline = 40;
    end
    
    methods
        function obj = TrimBaseline( prevJob )
           obj.name = 'Trim Pre/Post Baseline';
           if nargin > 0
               obj.prevJob = prevJob;
           end
        end
        
        function data = execute( obj, data )
            for i = 1:length(data)
                
                d = data(i).data;
                t = data(i).time;
                
                % check stims
                stims = data(i).stimulus.values;
                s = zeros(size(t,1),length(stims));
                for j = 1:length(stims)
                    s(:,j) = stims{j}.getStimVector( t );
                end
                
                s = abs( sum(s,2) );
                                
                t_min = t( find(s>0,1,'first') ) - obj.preBaseline;
                t_max = t( find(s>0,1,'last' ) ) + obj.postBaseline;
                
                lst = t >= t_min & t <= t_max;
                t = t(lst);
                d = d(lst,:);
                
                data(i).data = d;
                data(i).time = t;
                
            end
        end
        
        function options = getOptions( obj )
            options = [];
        end
           
        function obj = putOptions( obj, options )
        end
    end
    
end

