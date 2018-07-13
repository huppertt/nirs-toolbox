classdef Smooth < nirs.modules.AbstractModule
    %% Normalizes to the DC intensity
    %
    %
    properties
        sigma;
    end
    
    methods
        
        function obj = Smooth( prevJob )
            if nargin > 0, obj.prevJob = prevJob; end
            
            obj.name = 'Surface based smoothing';
            obj.sigma=6;
        end
        
        function data = runThis( obj, data )
           
            for i = 1:length(data)
               mesh=data(i).mesh.mesh;
               mesh.nodes=mesh.nodes(data(i).mesh.link.vertex,:);
               
               mesh=data(i).mesh.mesh;
               smoother=speye(size(mesh.nodes,1),size(mesh.nodes,1));
               lst=sub2ind(size(smoother),mesh.faces(:,1),mesh.faces(:,2));
               smoother(lst)=exp(-sum((mesh.nodes(mesh.faces(:,1),:)-mesh.nodes(mesh.faces(:,2),:)).^2,2)/obj.sigma^2);
               lst=sub2ind(size(smoother),mesh.faces(:,1),mesh.faces(:,3));
               smoother(lst)=exp(-sum((mesh.nodes(mesh.faces(:,1),:)-mesh.nodes(mesh.faces(:,3),:)).^2,2)/obj.sigma^2);
               lst=sub2ind(size(smoother),mesh.faces(:,2),mesh.faces(:,3));
               smoother(lst)=exp(-sum((mesh.nodes(mesh.faces(:,2),:)-mesh.nodes(mesh.faces(:,3),:)).^2,2)/obj.sigma^2);
               
               lst=sub2ind(size(smoother),mesh.faces(:,2),mesh.faces(:,1));
               smoother(lst)=exp(-sum((mesh.nodes(mesh.faces(:,2),:)-mesh.nodes(mesh.faces(:,1),:)).^2,2)/obj.sigma^2);
               lst=sub2ind(size(smoother),mesh.faces(:,3),mesh.faces(:,1));
               smoother(lst)=exp(-sum((mesh.nodes(mesh.faces(:,3),:)-mesh.nodes(mesh.faces(:,1),:)).^2,2)/obj.sigma^2);
               lst=sub2ind(size(smoother),mesh.faces(:,3),mesh.faces(:,2));
               smoother(lst)=exp(-sum((mesh.nodes(mesh.faces(:,3),:)-mesh.nodes(mesh.faces(:,2),:)).^2,2)/obj.sigma^2);
               
               
               smoother=smoother^6; 
               smoother=triu(smoother)+triu(smoother,1)';
               smoother=smoother./(sum(smoother,2)*ones(1,size(smoother,2)));
               smoother=smoother(data(i).mesh.link.vertex,data(i).mesh.link.vertex);
               d=data.data*data.projectors';
               d=d*smoother'; 
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