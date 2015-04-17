classdef ChannelStats
    %CHANNELSTATS Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        names           % variable names
        beta            % nconds x nchannels
        covb            % nconds x nconds x nchannels
        dfe          	% degrees of freedom
        demographics    % subject demographics
        probe           % probe geometry
    end
    
    methods
        function S = ttest(obj, C)
            % C -- matrix of contrast vectors
            % need to return beta, T, p, tcrit, names
            % default C to Identity
            beta = C*obj.beta;
            for i = 1:size(beta,2)
                covb(:,:,i) = C*obj.covb(:,:,i)*C';
                tstat(:,i)  = beta(:,i) ./ sqrt( diag(covb(:,:,i)) );
            end
            
            S.beta  = beta;
            S.covb  = covb;
            
            S.tstat = tstat;
            S.p     = tcdf(-abs(tstat), obj.dfe)/2;
            S.ppos  = 2*S.p .* (tstat >= 0);
            S.pneg  = 2*S.p .* (tstat <= 0);
            
            S.dfe   = obj.dfe;
            
            S.names = obj.transformNames(C);
        end
        
        function p = jointHypothesisTest(obj, h) %%% this wont work; monte carlo sampling from T-dist maybe?
            % M -- binary mask of F-test hypothesis
            % need to return F, p, df, fcrit
            error('')
        end
        
        function draw(obj, name, type, pthresh)
            obj.probe.draw()
        end
    end
    
    methods (Access = protected)
        function newNames = transformNames( obj, T )
            names = obj.names;
            for i = 1:size(T,1)
                newNames{i} = '';
                for j = 1:size(T,2)
                    c = T(i,j);
                    if c == 1
                        newNames{i} = [newNames{i} '+' names{j}];
                    elseif c == -1
                        newNames{i} = [newNames{i} '-' names{j}];
                    elseif c > 0
                        newNames{i} = [newNames{i} '+' num2str(c) names{j}];
                    elseif c < 0
                        newNames{i} = [newNames{i} num2str(c) names{j}];
                    end
                end
                
                if newNames{i}(1) == '+'
                    newNames{i} = newNames{i}(2:end);
                end
            end

        end
    end
    
end

