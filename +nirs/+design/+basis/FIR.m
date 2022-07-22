classdef FIR
    
    properties
        nbins = '16s';
        binwidth   = 1;
        isIRF = true;  
    end
    
    methods
        function [out,varargout] = convert( obj, s, t )
%             

%             nlag = round(Fs * obj.duration);
%             
%             out = lagmatrix(s, 0:nlag);
%             out(isnan(out)) = 0;
            
            
            if(~obj.isIRF)
                on = diff([0; s]) > 0;
                on=kron(on,ones(obj.binwidth,1));
                on(length(s)+1:end)=[];
            else
                on = s;
            end
            
            on = [on( floor(obj.binwidth/2)+1:end ); zeros(floor(obj.binwidth/2),1)];
            
            if(isstr(obj.nbins))
               Fs = 1/(t(2)-t(1));
                if(~isempty(strfind(obj.nbins,'s')))
                    nsec = str2num(obj.nbins(1:strfind(obj.nbins,'s')-1));
                    obj.nbins=ceil(nsec*Fs);
                else(~isempty(strfind(lower(obj.nbins),'fs*')))
                    nsec = str2num(obj.nbins(strfind(lower(obj.nbins),'fs*')+3:end));
                    obj.nbins=ceil(nsec*Fs);
                end
            end
            if(iscell(obj.nbins))
                for id=1:length(obj.nbins)
                    if(isstr(obj.nbins{id}))
                        Fs = 1/(t(2)-t(1));
                        if(~isempty(strfind(obj.nbins{id},'s')))
                            nsec = str2num(obj.nbins{id}(1:strfind(obj.nbins{id},'s')-1));
                            nbins(id)=ceil(nsec*Fs);
                        else(~isempty(strfind(lower(obj.nbins{id}),'fs*')))
                            nsec = str2num(obj.nbins{id}(strfind(lower(obj.nbins),'fs*')+3:end));
                            nbins(id)=ceil(nsec*Fs);
                        end
                    end
                end
                obj.nbins=nbins;
            end
           if(length(obj.nbins)>1 && obj.nbins(1)<0)
               on=[on; zeros(-obj.binwidth*obj.nbins(1),1)];
               n=sum(abs(obj.nbins))+1;
           else
               n=obj.nbins;
           end
            n=ceil(n);   
            
            f = kron(eye(n), ones(obj.binwidth,1));
            
            for i = 1:size(f,2)
               out(:,i) = filter(f(:,i), 1, on); 
            end
            
            out =[zeros(floor(obj.binwidth/2),size(out,2)); out(1:end-floor(obj.binwidth/2),:)];
      
           if(length(obj.nbins)>1 && obj.nbins(1)<0)
               out(1:-obj.binwidth*obj.nbins(1),:)=[];
           end
            
            
           if(nargout>1)
               varargout{1}=obj.nbins;
           end
           
        end
        
    end
    
end

