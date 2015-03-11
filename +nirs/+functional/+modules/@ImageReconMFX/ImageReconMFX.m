classdef ImageReconMFX < nirs.functional.AbstractModule
    %UNTITLED3 Summary of this class goes here
    %   Detailed explanation goes here
  
    properties
        formula = 'beta ~ cond*group + (1|subject)';
        jacobian = nirs.HashTable(); % key is subject name or "default"
        transMtx;
        mesh;
    end
    
    methods

        function obj = ImageReconMFX( prevJob )
           obj.name = 'Image Recon w/ Random Effects';
           if nargin > 0
               obj.prevJob = prevJob;
           end
           
           p = fileparts( which('nirs.functional.modules.ImageReconMFX') );
           load([p filesep 'wavelet_matrix.mat']);
           obj.transMtx = blkdiag(W,W);
        end
        
        function G = execute( obj, subj_stats )
            
            %% model table
            demo = nirs.functional.createDemographicsTable( subj_stats );
            
            [names, idx] = nirs.functional.getStimNames( subj_stats );
            tbl = [table(idx,names,'VariableNames',{'file','cond'}) ...
                demo(idx,:)];
            
            % parse table
            [x, z, names] = nirs.functional. ...
                parseWilkinsonFormula( obj.formula, tbl );
            
            x = sparse(x);
            z = sparse(z);
            
            %% Wavelet ( or other tranform )
            % the images are [left:hbo right:hbo left:hbr right:hbr]
            % there will be four diagonal blocks in W (there are already two)
            W = obj.transMtx;
            W = blkdiag(W,W);
            
            %% SVDS; PER FWD MODEL
            S = nirs.HashTable();
            U = nirs.HashTable();
            V = nirs.HashTable();
            for i = 1:length( obj.jacobian.keys )
                J = obj.jacobian.values{i};
                
                key = obj.jacobian.keys{i};
                [u,s,v] = svd(full([J.hbo J.hbr]),'econ');
                S( key ) = s;
                V( key ) = v;
                U( key ) = u;
                
            end
            clear J
            
            %% MASK; ACROSS FWD MODELS            
            mask = zeros( size(v,1), 1) > 0;
            for i = 1:length( obj.jacobian.keys )
                key = obj.jacobian.keys{i};
                
                % reparameterize model to fit wavelet coefs
                xiw = U(key)*S(key)*pinv(W*V(key));
                               
                % only fit coefs that have "enough" sensitivity
                lst = sqrt(sum(xiw.^2,1))';
                lst = abs(lst) > 1e-2*max(abs(lst));
                
                % we do the union across all forward models
                mask = mask | lst;
            end
            
            %% CALCULATE WAVELET MODEL
            XiW = nirs.HashTable();
            for i = 1:length( obj.jacobian.keys )
                key = obj.jacobian.keys{i};
                
                % reparameterize model to fit wavelet coefs
                xiw = U(key)*S(key)*pinv(W*V(key));
                               
                XiW(key) = xiw(:,mask);
            end
            clear xiw;
            
            %% ASSEMBLE X & Z
            X = []; Z = [];
            for i = 1:size(x,1)
                sname = demo.subject(1);
                
                if obj.jacobian.iskey( sname )
                    key = sname;
                elseif obj.jacobian.iskey('default')
                    key = 'default';
                else
                    error(['No forward model for subject: ' sname '.'])
                end
                
                X = [X; kron(x(i,:), XiW(key))];
                
                if ~isempty(z)
                    Z = [Z; kron(z(i,:), XiW(key))];
                end

            end
                        
            %% ORGANIZE DATA
            y = []; L = sparse([]);
            for i = 1:length(subj_stats)
                nCond = length( subj_stats(i).stimulus.keys );
                
                b = subj_stats(i).beta(1:nCond,:);
                covb = subj_stats(i).covb(1:nCond,1:nCond,:);
                
                l = sparse([]);
                for j = 1:size(covb,3)
                   l = blkdiag( l, chol(pinv(covb(:,:,j)),'lower') ); 
                end
                
                % permute matrix
                idx = [];
                for j = 1:nCond
                    idx = [idx  j:nCond:numel(b)];
                end
                l = l(idx,idx);

                % append
                b = b';
                y = [y; b(:)];
                L = blkdiag(L,l);
            end
            
            %% WHITEN DATA
            y = L*y; X = L*X; 
            if ~isempty( z )
                Z = L*Z;
                ZZ = Z*Z'; clear Z
            else
                ZZ = 0;
            end
            
            %% FITTING
            
            
           	vw = 1; vw0 = 1e16; % prior variance on wavelet coefs
            vz = 1; vz0 = 1e16; % prior variance on rfx
            
            X = full(X); ZZ = full(ZZ);
            
            iter = 0;
            while abs( (vw0-vw)/vw0 ) > 1e-2 && iter < 10
                
                vw0 = vw; vz0 = vz;

                Q = vz*ZZ + eye(size(X,1));
                
                iR = inv( chol(X*X'*vw + Q) );
                
                iX = vw*X'*(iR*iR');
                
                what = iX*y;

                H = X*iX;

                dfw = trace(H'*H);

                vw = (what'*what)/dfw;

                vz = ((y-X*what)'*(y-X*what) + trace(ZZ))/size(y,1);

                iter = iter+1;
                
                disp(iter);
            end
            clear Q H
            
            %% CONVERT BACK TO IMAGE SPACE
            ix = sparse( size(mask,1)*size(x,2),size(y,1) );
%             ix( repmat(mask,[size(x,2) 1]),: ) = iX;
            
            bhat = zeros(size(ix,1),1);
%             bhat( repmat(mask,[size(x,2) 1]) ) = what;
            
            nb = size(bhat,1); nw = size(W,1); ni = nb/nw;
            for i = 1:size(x,2)
                sname = demo.subject(1);
                
                if obj.jacobian.iskey( sname )
                    key = sname;
                elseif obj.jacobian.iskey('default')
                    key = 'default';
                else
                    error(['No forward model for subject: ' sname '.'])
                end
                
                iW = V(key)*S(key)'*U(key)'*XiW(key);
                
                idx1 = ((i-1)*nw+1):i*nw;
                idx2 = ((i-1)*sum(mask)+1):i*sum(mask);
                bhat( idx1 ) = iW*what( idx2 );
                
                ix(idx1,:) = iW*iX(idx2,:);
            end
            
%             for i = 2:size(x,2)
%                 iW = blkdiag(iW,iW);
%             end
%            
%             ix = iW*what;
% 
% %             ix = V*(pinv(W*V)*ix);
%             bhat = iX*y;


            %% PUT STATS
            lst = abs(bhat) < 0.01*max(abs(bhat));
            bhat(lst) = 0;
            
            G.iX        	= ix; clear ix;
            G.se            = sqrt( sum(G.iX.^2,2) );
            G.se( lst )     = Inf;
            G.tstat         = bhat./G.se;
            G.bhat          = bhat;
            G.dfe           = ceil(dfw);
            G.names         = names;
            
            
        end
        
        function options = getOptions( obj )
            options = [];
        end
           
        function obj = putOptions( obj, options )
        end
        
    end
    
end

