classdef ShortDistanceFilterImage < nirs.modules.AbstractModule
%% ShortDistanceFilter - This function filters based on short distance measurements
%
% Options:
%     ncomp - % number of components to remove

    properties
         maxnumcomp=6;
         sigma=60;
         braindepth=8;
         channelLst=[];
    end
    
    methods

        function obj =  ShortDistanceFilterImage( prevJob )
           obj.name = 'Short Distance Correction';         
           if nargin > 0
               obj.prevJob = prevJob;
           end
        end
        
        function data = runThis( obj, data )
            probe=data(1).probe;
            % resample data
            
            minX = min(probe.optodes.X);
            maxX = max(probe.optodes.X);
            dX = (maxX-minX)/10;
            minY = min(probe.optodes.Y);
            maxY = max(probe.optodes.Y);
            dY = (maxY-minY)/10;
            
            [X,Y,Z]=meshgrid([minX-dX*2:dX:maxX+dX*2],[minY-dY*2:dY:maxY+dY*2],[-obj.braindepth]);
            
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
            FwdModel.prop=nirs.media.tissues.brain(lambda,.7,50);
            FwdModel.Fm=0;
            
            FwdModel.probe=probe;
            Jacob=FwdModel.jacobian('spectral');
            L=[Jacob.hbo Jacob.hbr];
              
          n=length(X(:));
          Xn = repmat(X(:),1,n);
          Yn = repmat(Y(:),1,n);
          Zn = repmat(Z(:),1,n);
            
           Xn=Xn-Xn';
           Yn=Yn-Yn';
           Zn=Zn-Zn';
           Dist = Xn.^2 + Yn.^2 + Zn.^2;
            
           if(isempty(obj.channelLst))
               if(ismember('ShortSeperation',data(1).probe.link.Properties.VariableNames) && ...
                       any(data(1).probe.link.ShortSeperation))
                    channelLst=find(data(1).probe.link.ShortSeperation);
               else
                  channelLst=1:size(L,1);
               end
           elseif(isa(obj.channelLst,'table'))
               channelLst=find(ismember(probe.link,obj.channelLst));
           else
               channelLst=obj.channelLst;
           end
            
           smoother = exp(-Dist/obj.sigma^2);
           smoother=smoother/normest(smoother);
           
           obj.maxnumcomp=min(obj.maxnumcomp,length(channelLst));
           
           LSS=L(channelLst,:)*blkdiag(smoother,smoother);
           [U,S,V]=nirs.math.mysvd(LSS);
           LSS=U(:,1:obj.maxnumcomp)*S(1:obj.maxnumcomp,1:obj.maxnumcomp)*V(:,1:obj.maxnumcomp)';
           
           LSS=LSS/normest(LSS);
            
            
            
            for i = 1:numel(data)
                d = data(i).data;
                
                % remove mean
                m = mean(d,1);
                d = bsxfun(@minus, d, m);
                
                w=inv(chol(cov(d(:,channelLst))));
                wLSS=w*LSS;
                SStime=orth((L*blkdiag(smoother,smoother)*pinv(wLSS'*wLSS)*wLSS'*w*d(:,channelLst)')');
                SStime(:,end+1)=1;
                
                b=nirs.math.ar_irls(d,SStime,data(i).Fs*4);
                 
                d=d-SStime*b.beta;
                
                % add mean back
                d = bsxfun(@plus, d, m);
                
                % put back
                data(i).data = d;
                disp(['done ' num2str(i) ' of ' num2str(numel(data))]);
            end
        end
    end
    
end