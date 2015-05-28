classdef Anova < nirs.modules.AbstractModule
    %UNTITLED3 Summary of this class goes here
    %   Detailed explanation goes here
  
    properties
        formula = 'beta ~ 1 + group*cond + (1|subject)';
    end

    methods
        function obj = Anova( prevJob )
           obj.name = 'Anova Model';
           if nargin > 0
               obj.prevJob = prevJob;
           end
        end
        
        function G = runThis( obj, S )
            
            % demographics info
            demo = nirs.createDemographicsTable( S );
            
            % preallocate group stats
            G = nirs.core.AnovaStats();
            
            %% assemble table
            tbl = table();
            for i = 1:length(S)
                nCond = length(S(i).names);
                tbl = [tbl; [table(S(i).names(:),'VariableNames',{'cond'}) repmat(demo(i,:),[nCond 1])]];
            end
            
            %% loop through channels and fite mfx model
            for iChan = 1:size(S(1).probe.link.source,1)
                
                % get hemodynamic response and covariance
                beta = []; W = sparse([]);
                for i = 1:length(S)
                    % coefficients
                    beta = [beta; S(i).beta(1:nCond,iChan)];
                    
                    % design whitening transform from svd
                    [u, s, ~] = svd( S(i).covb(1:nCond, 1:nCond, iChan) );
                    s = 1./diag(sqrt(s));
                    w = diag(s)*u';
                    
                    % put them in giant block diag matrix
                    W = blkdiag(W, w);
                end
                
                % weighted fit to get design matrices
                lm1 = fitlme([table(beta) tbl], obj.formula, ...
                    'FitMethod', 'REML', 'Weights', full(sqrt(diag(W'*W))));
                                
                a = lm1.anova();
                
                G.F(:,iChan) = a.FStat;
                
            end
            
            G.df1   = a.DF1;
            G.df2   = a.DF2;
            G.names = a.Term;
            G.probe = S(1).probe;
            
        end
    end
    
end

