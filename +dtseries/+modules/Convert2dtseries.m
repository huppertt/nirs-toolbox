classdef Convert2dtseries < nirs.modules.AbstractModule
    % This function converts NIRS, EEG, or MEG (WIP) data to a dense time
    % series variable
    
    properties
        jacobian = Dictionary(); % key is subject name or "default"
        basis;  % Basis set from nirs.inverse.basis
        mask = [];   % Mask (e.g. cortical contraint)
        prior = Dictionary();  % The prior on the image; default = 0 (Min Norm Estimate)
        mesh;  % Reconstruction mesh (needed to create the ImageStats Class)
        probe = Dictionary();  % This is the probe for the Jacobian model.  These needs to match the jacobian
    end
    
    methods
        
        function obj = Convert2dtseries( prevJob )
            obj.name = 'Image Recon to dense time series';
            if nargin > 0
                obj.prevJob = prevJob;
            end
            
            %             nVox=20484;
            %             obj.basis=nirs.inverse.basis.identity(nVox);
            %
            %             prior.hbo=zeros(nVox,1);
            %             prior.hbr=zeros(nVox,1);
            %             obj.prior('default')=prior;
            
        end
        
        function G = runThis( obj,data )
            
            
            % Mask of the reconstruction volume
            if(~isempty(obj.mask))
                LstInMask=find(obj.mask);
            else
                LstInMask=1:size(obj.basis.fwd,1);
            end
            
            %% Wavelet ( or other tranform )
            Basis = obj.basis.fwd; % W = W(1:2562,:);
            
            %Let's make the forward models
            L = obj.jacobian;
            Lfwdmodels=Dictionary();
            Probes=Dictionary();
            
            for i = 1:L.count
                key = L.keys{i};
                J =L(key);
                flds=fields(J);
                
                if(~ismember(key,obj.probe.keys))
                    Probes(key)=obj.probe('default');
                else
                    Probes(key)=obj.probe(key);
                end
                l=[];
                for j=1:length(flds)
                    l=setfield(l,flds{j},(J.(flds{j})(:,LstInMask).*(J.(flds{j})(:,LstInMask)>10*eps(1)))*...
                        Basis(LstInMask,:));
                end
                Lfwdmodels(key)=l;
            end
            
            % Do a higher-order generalized SVD
            [US,V]=nirs.math.hogSVD(Lfwdmodels.values);
            % [U1,U2...,V,S1,S2...]=gsvd(L1,L2,...);
            % L1 = U1*S1*V'
            % L2 = U2*S2*V'
            
            % Store back into the forward model
            Lfwdmodels.values=US;
            
            
            %Make sure the probe and data link match
            for idx=1:length(data)
                sname=data(idx).demographics('subject');
                if(~isempty(sname) && ismember(sname,Lfwdmodels.keys))
                    key=sname;
                else
                    key='default';
                end
                thisprobe=Probes(key);
                [ia,ib]=ismember(thisprobe.link,data(idx).probe.link,'rows');
                
                data(idx).data=data(idx).data(:,ib);
                data(idx).probe.link=data(idx).probe.link(ib,:);
                
            end
            
            
            
            d=[];
            W=zeros(size(data,2));
            for idx=1:length(data)
                
                 W=W+cov(data(i).data);
                
                key = data(i).demographics('subject');
                if(~Lfwdmodels.iskey(key))
                    key='default';
                end
                xx=[];
                
                xlocal=[];
                for fIdx=1:length(flds)
                    L=Lfwdmodels(key);
                    x2=L.(flds{fIdx});
                    
                    xlocal=[xlocal x2];
                end
               d=[d data(idx).data'];
                %X=[X xlocal];
            end
            [q,r]=qr(d,0);
            X=q'*xlocal;
            d=r;
            
            [u, s, ~] = svd(W, 'econ');
            W=diag(1./diag(eps(1)+sqrt(s))) * u';
            
            d=W*d;
            X=W*X;
            for idx=1:size(X,2)
                x=X(:,idx);
                VarMDU(idx)=inv(x'*x+eps(1));
            end
            lst=find(any(reshape(sqrt(VarMDU),[],length(flds))<1/sqrt(eps(1))/10,2));
            
            n=length(lst);
            lst2=[];
            for id=1:length(flds)
                q=sparse(n*length(flds),n*length(flds));
                q((id-1)*n+1:id*n,(id-1)*n+1:id*n)=sparse(diag(ones(n,1)));
                Q{id}=q;
                lst2=[lst2; lst+(id-1)*length(VarMDU)/length(flds)];
            end

            
            
            
            R={speye(size(d,1),size(d,1))};
            
            Beta0= zeros(size(X,2),1);
            [lambda,Beta,Stats]=nirs.math.REML(d,X(:,lst2),Beta0(lst2),R,Q);
            
            V=obj.basis.fwd*V(:,lst);
            
           a=sum(abs(V),2);
           lst=zscore(a)<1;
           V(lst,:)=0;
           
            vertex=repmat([1:size(V,2)]',length(flds),1);
            type=reshape(repmat(flds',size(V,2),1),[],1);
            link=table(vertex,type);
            cnt=0;
            for idx=1:length(data)
                G(idx)=dtseries.core.Data;
                G(idx).description=data(idx).description;
                G(idx).demographics=data(idx).demographics;
                G(idx).stimulus=data(idx).stimulus;
                G(idx).time=data(idx).time;
                G(idx).data = Beta(:,cnt+[1:length(G(idx).time)])';
                G(idx).mesh=dtseries.core.Mesh(obj.mesh,link);
                G(idx).cov=Stats.tstat.covb;
                G(idx).projectors=V;
                cnt=cnt+length(G(idx).time);
            end
            
            
        end
    end
end