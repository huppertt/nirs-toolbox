classdef IterImageReconMFX < nirs.functional.AbstractModule
    %UNTITLED3 Summary of this class goes here
    %   Detailed explanation goes here
  
    properties
        formula = 'beta ~ cond*group + (1|subject)';
        jacobian = nirs.HashTable(); % key is subject name or "default"
        dummyCoding = 'full';
        transMtx;
        mesh;
    end
    
    methods

        function obj = IterImageReconMFX( prevJob )
           obj.name = 'Iterative Image Recon w/ Random Effects';
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
                parseWilkinsonFormula( obj.formula, tbl, 'true', obj.dummyCoding );
            
%             x = sparse(x);
%             z = sparse(z);
            
            %% Wavelet ( or other tranform )
            % the images are [left:hbo right:hbo left:hbr right:hbr]
            % there will be four diagonal blocks in W (there are already two)
            W = obj.transMtx; % W = W(1:2562,:);
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
            clear J u s v
            
            %% MASK; ACROSS FWD MODELS            
            % maskw = zeros( size(v,1), 1) > 0;
            maskb = zeros( size(W,1), 1) > 0;
            for i = 1:length( obj.jacobian.keys )
                key = obj.jacobian.keys{i};
                
                % reparameterize model to fit wavelet coefs
                X   = U(key)*S(key)*V(key)';
%                 xiw = U(key)*S(key)*pinv(W*V(key));
                               
%                 % only fit coefs that have "enough" sensitivity
%                 lst = sqrt(sum(xiw.^2,1))';
%                 lst = abs(lst) > 1e-2*max(abs(lst));
                
%                 % we do the union across all forward models
%                 maskw = maskw | lst;
                
                % only fit coefs that have "enough" sensitivity
                lst = sqrt(sum(X.^2,1))';
                lst = abs(lst) > 1e-2*max(abs(lst));
                
                % we do the union across all forward models
                maskb = maskb | lst;
            end
            
            %% ORGANIZE DATA
            y = {}; P = {};
            for i = 1:length(subj_stats)
                nCond = length( subj_stats(i).stimulus.keys );
                
                b = subj_stats(i).beta(1:nCond,:);
                covb = subj_stats(i).covb(1:nCond,1:nCond,:);
                
                l = sparse([]);
                for j = 1:size(covb,3)
                   l = blkdiag( l, inv(chol(covb(:,:,j)))); 
                end
                
                % permute matrix
                idx = [];
                for j = 1:nCond
                    idx = [idx  j:nCond:numel(b)];
                end
                l = l(idx,idx);

                % append
                y{i} = b;
                P{i} = l*l';
            end
            
            vec = @(x) x(:);
            
            %% ITERATION
            m = zeros(size(W,1),size(x,2)); % group mean
            v = 100*ones( size(m) );        % group variance
            for iter = 1:10
                m0 = m; v0 = v;

                %% FIT EACH FILE
                for iFile = 1:size(y)
                    idx = x(iFile,:)>0;
                    b0  = vec( m(:,idx) );
                    vb0 = vec( v(:,idx) );
                    
                    
                end
                
                %% AVERAGE ACROSS FILES
                
                
                
            end
           
            
            %% PUT STATS
            lst = abs(bhat) < 0.01*max(abs(bhat))  | repmat(~maskb,[length(bhat)/length(maskb) 1]);
            bhat(lst) = 0;
            
            G.iX        	= ix; clear ix;
            G.se            = sqrt( sum(G.iX.^2,2) );
            G.se( lst )     = Inf;
            G.tstat         = bhat./G.se;
            G.bhat          = bhat;
            G.dfe           = ceil(dfw);
            G.names         = names;
            G.p             = 2*tcdf(-abs(G.tstat),G.dfe);
            G.ppos        	= tcdf(-G.tstat,G.dfe);
            G.pneg        	= tcdf(G.tstat,G.dfe);
            
            
        end
        
        function options = getOptions( obj )
            options = [];
        end
           
        function obj = putOptions( obj, options )
        end
        
    end
    
end

