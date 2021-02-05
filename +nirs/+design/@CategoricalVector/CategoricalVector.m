classdef CategoricalVector
    %This class holds categorical vectors
    
    properties
        name
        values;
        time
    end
    
    methods
        function vec = getStimVector( obj, time )
            
            [~,i]=ismember(obj.values,categories(obj.values));
            
            vec = interp1( obj.time, i, time ,'nearest','extrap');
            c=categories(obj.values);
            vec=categorical([c(vec)]);
          
            
        end
        
        function draw(obj)

            figure;
            plot(obj.time,obj.values);
            xlabel('Time (sec)')
            ylabel(obj.name);

            
        end
        
    end
    
end