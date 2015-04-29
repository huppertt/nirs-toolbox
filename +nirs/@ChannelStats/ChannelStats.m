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
            
            if nargin == 1
                C = eye(size(obj.beta,1));
            end
            
            beta = C*obj.beta;
            for i = 1:size(beta,2)
                covb(:,:,i) = C*obj.covb(:,:,i)*C';
                tstat(:,i)  = beta(:,i) ./ sqrt( diag(covb(:,:,i)) );
            end
            
            S.beta  = beta;
            S.covb  = covb;
            
            S.tstat = tstat;
            S.p     = tcdf(-abs(tstat), obj.dfe) * 2;
            S.ppos  = S.p .* (tstat >= 0) / 2;
            S.pneg  = S.p .* (tstat <= 0) / 2;
            
            S.dfe   = obj.dfe;
            
            S.names = obj.transformNames(C);
        end
        
        function S = ftest(obj, m)
            % this uses Hotelling's T squared test for joint 
            % hypothesis testing
            
            if nargin == 1
                m = ones(size(obj.beta,1),1) > 0;
            end
            
            if ~islogical(m)
                m = m > 0;
                warning('Converting mask to true/false.')
            end
                        
            n = obj.dfe;
            k = size(obj.beta,1);
            
            for i = 1:size(obj.beta,2)
                b = m(:) .* obj.beta(:,i);
                T2(i,1)     = b'*pinv(obj.covb(:,:,i))*b;
                F(i,1)      = (n-k) / k / (n-1) * T2(i);
                p(i,1)      = fcdf(1/F(i), n-k, k);
            end
            
            S.T2    = T2;
            S.F     = F;
            S.p     = p;
            S.df2   = n-k;
            S.df1   = k;
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

