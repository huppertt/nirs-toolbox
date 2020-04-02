classdef Normalize < nirs.modules.AbstractModule
    %% Normalizes to the DC intensity
    %
    %
    properties
        
    end
    
    methods
        
        function obj = Normalize( prevJob )
            if nargin > 0, obj.prevJob = prevJob; end
            
            obj.name = 'Dense time series DC normalization';
            
        end
        
        function data = runThis( obj, data )
           
            for i = 1:length(data)
               
                d = data(i).data*data(i).projectors';
                d = d./(ones(size(d,1),1)*mean(d,1))-1;
                
                lst=find([any(isnan(d),1) | all(d==0,1) | sqrt(var(d,[],1))<eps(1)*10]);
                d(:,lst)=[];
                data(i).projectors(lst,:)=[];
                data(i).mesh.link(lst,:)=[];
                
                if(size(d,2)>size(d,1))
                    [u,s,proj]=nirs.math.mysvd(d);
                     lst=find(diag(s)<eps(single(1)));
                     u(:,lst)=[]; s(lst,:)=[]; s(:,lst)=[]; proj(:,lst)=[];
                    data(i).data=u*s;
                     a=data.projectors'*proj;
                    data(i).cov = a'*data.cov*a;
                    
                    data(i).projectors=sparse(proj);
                else
                    data(i).data=d*data(i).projectors;
                end
                
               
            end
            
        end
        
        
        
    end
    
end

