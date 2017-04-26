classdef General
% General HRF response
    
    properties
        IRF;
        time;
    end
    
    methods
        function obj=General()
            fs=20;
            b=nirs.design.basis.Canonical;
            obj.time=[0:1/fs:40]';
            obj.IRF=b.convert([1; zeros(length(obj.time)-1,1)],obj.time);
        end
        function out = convert( obj, s, t )
        
          
            h=interp1(obj.time,obj.IRF,t,'linear');
            h(isnan(h))=0;
           
            % convert stim vectors
            for i=1:size(h,2)
                out(:,i) = filter(h(:,i), 1, s);
            end
     
        end
    end
    
            
end