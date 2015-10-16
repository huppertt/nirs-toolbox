classdef FIR
    
    properties
        nbins = 10;
        binwidth   = 5;
        isIRF = true;  
    end
    
    methods
        function out = convert( obj, s, t )
%             
%             Fs = 1/(t(2)-t(1));
%             nlag = round(Fs * obj.duration);
%             
%             out = lagmatrix(s, 0:nlag);
%             out(isnan(out)) = 0;

            if(~obj.isIRF)
                on = diff([0; s]) > 0; 
            else
                on = s;
            end
            on = [on( floor(obj.binwidth/2)+1:end ); zeros(floor(obj.binwidth/2),1)];
            
            f = kron(eye(obj.nbins), ones(obj.binwidth,1));
            
            for i = 1:size(f,2)
               out(:,i) = filter(f(:,i), 1, on); 
            end
            
             out =[zeros(floor(obj.binwidth/2),size(out,2)); out(1:end-floor(obj.binwidth/2),:)];
      
        end
        
    end
    
end

