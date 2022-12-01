classdef FIR_smoothing
    
    properties
        nbins = '16s';
        binFWHM = '4s';
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
                on(length(s)+1:end)=[];
            else
                on = s;
            end
            
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
            
            if(isstr(obj.binFWHM))
                FWHM = str2num(obj.binFWHM(1:strfind(obj.binFWHM,'s')-1));
            else
                FWHM=obj.binFWHM;
            end
            
            tt=[-FWHM*8:1/Fs:FWHM*8]';
            sig=FWHM/(2*sqrt(2*log(2)));
            G = 1/(sig*sqrt(2*pi))*exp(-tt.^2/sig^2/2)/Fs;
            on=[on; zeros(length(G),1)];
            
            for i = 1:n
               out(:,i) = filter([zeros(i-1,1); G], 1, on); 
            end
             out(find(tt<0),:)=[];
             out=out(1:length(t),:);
            
           if(nargout>1)
               varargout{1}=obj.nbins;
           end
           
        end
        
    end
    
end
