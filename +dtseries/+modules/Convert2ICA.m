classdef Convert2ICA < nirs.modules.AbstractModule
    %% Convert2ICA - Converts a dense time series to a ICA projector model.
    %
    %
    properties
        verbose;
    end
    
    methods
        
        function obj = Convert2ICA( prevJob )
            if nargin > 0, obj.prevJob = prevJob; end
            
            obj.name = 'Dense time series ICA projection';
            obj.verbose=false;
        end
        
        function data = runThis( obj, data )
            vec = @(x) x(:);
            
            for i = 1:length(data)
                
                ii=unique(data(i).mesh.link.type);
                mixedsig=[];
                dd=[];
                for idx=1:length(ii)
                    mixedsig=[mixedsig; data(i).data(:,ismember(data(i).mesh.link.type,ii(idx)))*data(i).projectors'];
                    dd=[dd; data(i).data(:,ismember(data(i).mesh.link.type,ii(idx)))];
                end
                
                dd=dd-ones(size(dd,1),1)*mean(dd,1);
                [u,s,v]=nirs.math.mysvd(cov(dd));
                [u2,s2,v2]=nirs.math.mysvd(data(i).projectors);
                
                dwm=u2*s2*v2'*u*sqrt(s);
                %u'*v2*u2'*dwm=s2*sqrt(s);
                wm=diag(1./sqrt(diag(s)))*u'*v2*diag(1./diag(s2))*u2';
                nv = wm*mixedsig';
                
                if(~obj.verbose)
                    verbose='off';
                else
                    verbose='on';
                end
                [A, W] = fpica(nv, wm, dwm,'symm',size(dd,2),...
                    'tanh','off',1,1,1,'on',0.001,5000,100,...
                    'rand',1,1,'off',1,verbose);
                
                P=data(i).projectors'*W';
                for idx=1:length(ii)
                    lst=ismember(data(i).mesh.link.type,ii(idx));
                    data(i).data(:,lst)= data(i).data(:,lst)*P;
                end
                data(i).projectors=A;
                
                % print progress
                obj.printProgress( i, length(data) )
            end
            
        end
        
        
        
    end
    
end

