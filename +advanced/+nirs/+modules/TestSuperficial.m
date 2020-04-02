classdef TestSuperficial < nirs.modules.AbstractModule
    %% TestSuperficial - preforms "ROI" based testing of each channel against the superficial layer.
    %
    %
    properties
        braindepth = 10;  % depth of brain layer (in mm)
        sigma = 120;  % Spatial smoothing for superficial layer
        method = 2;
        %allow_flip=false;  % This prevents sign flips when the mean over all channels is non-zero.
    end
    
    methods
        function obj = TestSuperficial( prevJob )
            obj.name = 'Test Superficial';
            if nargin > 0
                obj.prevJob = prevJob;
            end
        end
        
        function data = runThis( obj, data )
            
            for i = 1:numel(data);
                data(i)=sorted(data(i),{'source','detector','type'});
            end
            
            % Assume all the probes are the same
            probe=data(1).probe;
            
            % Compute the optical forward model based on the slab model
            minX = min(probe.optodes.X);
            maxX = max(probe.optodes.X);
            dX = (maxX-minX)/10;
            minY = min(probe.optodes.Y);
            maxY = max(probe.optodes.Y);
            dY = (maxY-minY)/10;
            
            [X,Y,Z]=meshgrid([minX-dX:dX:maxX+dX],[minY-dY:dY:maxY+dY],[-obj.braindepth -obj.braindepth*2]);
            
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
            brainmask = (Z<-obj.braindepth);
            
            for i = 1:numel(data)
                probe=data(i).probe;
                c=eye(size(L,1));
                for idx=1:height(probe.link)
                    basisX = X-(probe.srcPos(probe.link.source(idx),1)+probe.detPos(probe.link.detector(idx),1))/2;
                    basisY = Y-(probe.srcPos(probe.link.source(idx),2)+probe.detPos(probe.link.detector(idx),2))/2;
                    basis = exp(-basisX.^2/obj.sigma^2).*exp(-basisY.^2/obj.sigma^2);
                    
                    lst=find(ismember(probe.link.type,probe.link.type(idx)));
                    
                   % w=inv(chol(data(i).covb(lst,lst)));
                         
                    [u, s, ~] = svd(data(i).covb(lst,lst), 'econ');
                    %W = blkdiag(W, diag(1./diag(sqrt(s))) * u');
                    w = pinv(s).^.5 * u';
                
                
                
                    if(obj.method==2)
                        % In this option, we are doing a "reconstruction" of the inverse problem using a
                        % smoothed/skin-resstricted basis set from the leave-one-out set of data.  Then, we do a T-test on the
                        % residual of the channel we left out.  E.g. is
                        % the channel of interest statstically different
                        % then that predicted from a skin response of all
                        % the other channels
                        
                        lambda=1./abs(data(i).tstat(idx))^2;
                        ll=w*L(lst,:)*[skinmask(:).*basis(:)];
                        ll(ismember(lst,idx),:)=0;
                        c(idx,lst) = c(idx,lst)*w-L(idx,:)*(skinmask(:).*basis(:))*pinv(ll'*ll+lambda*eye(size(ll,2)))*ll'*w;
                        
                        
                        %normalize to make sure it sums to zero
                        lstP=find(c(idx,:)>0);
                        lstN=find(c(idx,:)<0);
                        
                        c(idx,lstP)=3*c(idx,lstP)/sum(c(idx,lstP));
                        c(idx,lstN)=-c(idx,lstN)/sum(c(idx,lstN));
                    else
                        % In this example, we do an ROI based model and do
                        % a T-test of an ROI defined by a low spatial
                        % resolution skin basis and the brain region seen by the single detector
                        a=(L(lst,:)*(skinmask(:).*basis(:)))';
                        c(idx,lst) = c(idx,lst)-a;
                        
                        %normalize to make sure it sums to zero
                        lstP=find(c(idx,:)>0);
                        lstN=find(c(idx,:)<0);
                        
                        c(idx,lstP)=3*c(idx,lstP)/sum(c(idx,lstP));
                        c(idx,lstN)=-c(idx,lstN)/sum(c(idx,lstN));
                    end
                    
                end
                
                
                cc = kron(c,eye(length(data(i).conditions)));
                
                
                if(0) %obj.allow_flip)
                    % I don't think it ever makes sense to allow this
                    data(i).beta = cc * data(i).beta;
                    
                else
                    data(i).beta = (cc * data(i).beta).*(sign(data(i).beta)==sign(cc * data(i).beta));
                
                end
                s=min(diag(data(i).covb))/10;
                data(i).covb = cc * data(i).covb * cc';
                data(i).covb=data(i).covb+eye(size(cc,1))*(s);
                
            end
            
        end
    end
end
