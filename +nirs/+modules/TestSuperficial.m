classdef TestSuperficial < nirs.modules.AbstractModule
%% TestSuperficial - preforms "ROI" based testing of each channel against the superficial layer.
% 
%    
    properties
       braindepth = 10;  % depth of brain layer (in mm)
       sigma = 100;  % Spatial smoothing for superficial layer
       allow_flip=true;  % This prevents sign flips when the mean over all channels is non-zero.
    end
    
    methods
        function obj = TestSuperficial( prevJob )
           obj.name = 'Test Superficial';
           if nargin > 0
               obj.prevJob = prevJob;
           end
        end
        
        function data = runThis( obj, data )
        
                for i=1:length(data);
                    data(i)=sorted(data(i),{'source','detector','type'});
                end
            
                % Assume all the probes are the same
                probe=data(1).probe;
                
                % Compute the optical forward model based on the slab model
                minX = min(probe.optodes.X);
                maxX = max(probe.optodes.X);
                dX = (maxX-minX)/3;
                minY = min(probe.optodes.Y);
                maxY = max(probe.optodes.Y);
                dY = (maxY-minY)/3;
                
                [X,Y,Z]=meshgrid([minX-dX:dX/30:maxX+dX],[minY-dY:dY/30:maxY+dY],[-1:-2:-obj.braindepth*2]);
                
                mesh=nirs.core.Mesh;
                mesh.nodes=[X(:) Y(:) Z(:)];
                
                lambda=unique(probe.link.type);
                if(iscell(lambda));
                    disp('using 808nm as appromation to hemoglobin fwd-model')
                    lambda=808;
                    probe.link.type=repmat(lambda,height(probe.link),1);
                end;
                
                FwdModel=nirs.forward.ApproxSlab;
                FwdModel.mesh=mesh;
                FwdModel.prop=nirs.media.tissues.brain(.7,50,lambda);
                FwdModel.Fm=0;
                
                FwdModel.probe=probe;
                Jacob=FwdModel.jacobian;
                L=Jacob.mua; 
                L=L./(sum(L,2)*ones(1,size(L,2))); 
                skinmask = (Z>=-obj.braindepth);
                
                probe=data(1).probe;
                 c=eye(size(L,1));
                 for idx=1:height(probe.link)
                     basisX = X-(probe.srcPos(probe.link.source(idx),1)+probe.detPos(probe.link.detector(idx),1))/2;
                     basisY = Y-(probe.srcPos(probe.link.source(idx),2)+probe.detPos(probe.link.detector(idx),2))/2;
                     basis = exp(-basisX.^2/obj.sigma^2).*exp(-basisY.^2/obj.sigma^2);
                     
                     lst=find(ismember(probe.link.type,probe.link.type(idx)));
                     c(lst,idx) = c(lst,idx)-L(lst,:)*(skinmask(:).*basis(:));
                    
                     % normalize to make sure it sums to zero
                     lstP=find(c(:,idx)>0); 
                     lstN=find(c(:,idx)<0); 
                     
                     c(lstP,idx)=c(lstP,idx)/sum(c(lstP,idx));
                     c(lstN,idx)=-c(lstN,idx)/sum(c(lstN,idx));
                     
                 end

                for i=1:length(data)
                    cc = kron(c,eye(length(data(i).conditions)));
                    
                    if(obj.allow_flip)
                        data(i).beta = cc * data(i).beta;
                    else
                        data(i).beta = max(cc * abs(data(i).beta),0).*sign(data(i).beta);
                    end
                    
                    data(i).covb = cc * data(i).covb * cc'+eye(size(cc,1))*1E-12;

                end
                
        end
    end
end
