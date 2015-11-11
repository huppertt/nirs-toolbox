classdef ApproxSlab
    %SLABFORWARDMODEL Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        probe;
        prop;
        mesh;
        Fm = 110;
    end
    
    methods
        %% Constructor
        function obj = ApproxSlab( mesh, prop, probe, Fm )
            if nargin > 0, obj.mesh = mesh; end
            if nargin > 1, obj.prop = prop; end
            if nargin > 2, obj.probe = probe; end
            if nargin > 3, obj.Fm = Fm; end
        end
        
        %% Methods
        function meas = measurement( obj )
            
            dist=obj.probe.distances;
            m=zeros(size(dist));
            for iLambda = 1:length(obj.prop.lambda)
                V = obj.prop.v ./ obj.prop.ri;
                D = V ./ (3 * (obj.prop.musp(iLambda) +obj.prop.mua(iLambda)));
                K = sqrt(V .* obj.prop.mua(iLambda) ./ D - i * (2*pi * (obj.Fm * 1e6)) ./ D);
                
                % 2pt Green's function for photon density (freqency-domain)
                %phiS = exp(-K * src_r) ./ (4 * pi * D * src_r);
                lst=find(ismember(obj.probe.link.type,obj.prop.lambda(iLambda)));
                m(lst) = exp(-K * dist(lst)) ./ ( dist(lst));
            end
            meas = nirs.core.Data(m,...
                obj.probe,0,obj.Fm);
        end
        
        function [J,meas] = jacobian( obj,type )
            
            if nargin < 2
                isSpectral = false;
            elseif strcmpi( type,'standard' )
                isSpectral = false;
            elseif strcmpi( type,'spectral' );
                isSpectral = true;
            else
                error('Jacobian can either be ''standard'' or ''spectral''.')
            end
            
            types = unique( obj.probe.link.type );
            assert( isnumeric( types ) );
            [~,~,iType] = unique( obj.probe.link.type );
            
            meas=obj.measurement;
            mesh = obj.combinemesh;
            for idx=1:size(obj.probe.srcPos,1)
                src_r(idx,:)=sqrt((mesh.nodes(:,1)-obj.probe.srcPos(idx,1)).^2+...
                    (mesh.nodes(:,2)-obj.probe.srcPos(idx,2)).^2+...
                    (mesh.nodes(:,3)-obj.probe.srcPos(idx,3)).^2);
            end
            for idx=1:size(obj.probe.detPos,1)
                det_r(idx,:)=sqrt((mesh.nodes(:,1)-obj.probe.detPos(idx,1)).^2+...
                    (mesh.nodes(:,2)-obj.probe.detPos(idx,2)).^2+...
                    (mesh.nodes(:,3)-obj.probe.detPos(idx,3)).^2);
            end
            
            Jmua = zeros( size(obj.probe.link,1), size(mesh.nodes,1) );
            for iLambda = 1:length(obj.prop.lambda)
                
                V = obj.prop.v ./ obj.prop.ri;
                D = V ./ (3 * (obj.prop.musp(iLambda) +obj.prop.mua(iLambda)));
                K = sqrt(V .* obj.prop.mua(iLambda) ./ D - i * (2*pi * (obj.Fm * 1e6)) ./ D);
                
                % 2pt Green's function for photon density (freqency-domain)
                %phiS = exp(-K * src_r) ./ (4 * pi * D * src_r);
                phiS = exp(-K * src_r) ./ ( src_r);
                phiS(find(src_r==0))=1;
                
                phiD = exp(-K * det_r) ./ ( det_r);
                phiD(find(det_r==0))=1;
                lst=find(ismember(obj.probe.link.type,obj.prop.lambda(iLambda)));
                Jmua(lst,:)=phiS(obj.probe.link.source(lst),:).*phiD(obj.probe.link.detector(lst),:);
                
            end
            
            
            if ~isSpectral
                J.mua = Jmua;
            else
                % convert jacobian to conc
                ext = nirs.media.getspectra( types );
                
                ehbo = ext(iType,1);
                ehbr = ext(iType,2);
                
                J.hbo = bsxfun(@times,ehbo,Jmua);
                J.hbr = bsxfun(@times,ehbr,Jmua);
                
            end
            
        end
        
        
        function mesh = combinemesh(obj)
            mesh = nirs.core.Mesh;
            
            mesh=obj.mesh(1);
            
            for idx=2:length(obj.mesh)
                n=size(mesh.nodes,1);
                mesh.nodes=[mesh.nodes; obj.mesh(idx).nodes];
                mesh.faces=[mesh.faces; obj.mesh(idx).faces+n];
                mesh.elems=[mesh.elems; obj.mesh(idx).elems+n];
                mesh.regions=[mesh.regions; obj.mesh(idx).regions];
                try; mesh.fiducials=[mesh.fiducials; obj.mesh(idx).fiducials]; end;
            end
            
            
        end
    end
    
end
