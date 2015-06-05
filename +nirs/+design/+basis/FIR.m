classdef FIR
    
    properties
        npoints = 6;
        width   = 3;
    end
    
    methods
        function out = convert( obj, s, t )
            
%             Fs = 1/(t(2)-t(1));
%             nlag = round(Fs * obj.duration);
%             
%             out = lagmatrix(s, 0:nlag);
%             out(isnan(out)) = 0;

            on = diff([0; s]) > 0;
            
            f = kron(eye(obj.npoints), ones(obj.width,1));
            
            for i = 1:size(f,2)
               out(:,i) = filter(f(:,i), 1, on); 
            end
      
      
        end
        
    end
    
end

