classdef FIR
    
    properties
        npoints = 10;
        width   = 5;
    end
    
    methods
        function out = convert( obj, s, t )
            
%             Fs = 1/(t(2)-t(1));
%             nlag = round(Fs * obj.duration);
%             
%             out = lagmatrix(s, 0:nlag);
%             out(isnan(out)) = 0;

            on = diff([0; s]) > 0; 
            
            on = [on( floor(obj.width/2)+1:end ); zeros(floor(obj.width/2),1)];
            
            f = kron(eye(obj.npoints), ones(obj.width,1));
            
            for i = 1:size(f,2)
               out(:,i) = filter(f(:,i), 1, on); 
            end
      
      
        end
        
    end
    
end

