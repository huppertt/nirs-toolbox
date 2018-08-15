classdef ImageStats
    %% IMAGESTATS - Holds regression stats for channel space.
    % 
    % Properties: 
    %     description  - description of data (e.g. filename)
    %     mesh         - The reconstruction mesh
    %     variables    - a table containing w/ source, detector, type, cond for
    %                    each regression coefficient
    %     beta         - regression coefficients
    %     covb         - covariance of beta
    %     dfe          - degrees of freedom
    %     probe        - Probe object describing measurement geometry
    %     demographics - Dictionary containing demographics info
    %                    (e.g. demographics('age') returns 28)
    %     typeII_StdE  - StdErr used to compute the typeII 
    %     conditions   - (dependent) list of stim conditions
    %     tstat        - (dependent) t-stats of beta
    %     p            - (dependent) p-values of beta
    %     q            - (dependent) q-values of beta (false discovery rate)
    %       
    %  Methods:
    %     getCritT    - returns critical t value
    %     draw        - draws beta or tstat values on top of probe geometry
    %     table       - returns a table of all stats (minus full covariance)
    %     ttest       - performs a t-test and returns a new ChannelStats object
    %     ftest       - performs an F-test across conditions 
    %     jointTest   - performs a joint hypothesis test across all channels in
    %                   each SD pair (i.e. joint test of hbo & hbr for S1-D1)
    %     sorted      - returns new ChannelStats object with sorted channels

    properties
        description     % description of data (e.g. filename)
        mesh            % the reconstuction mesh
        variables       % table describing regression coefficients
        
        beta            % regression coefficients
        covb_chol       % covariance of beta =covb_chol*covb_chol'  (since storing the while thing is too large)
        dfe          	% degrees of freedom
        typeII_StdE     % StdErr used to compute the typeII 
        probe           % Probe object describing measurement geometry
        demographics    % Dictionary containing demographics info
    end
    
    properties ( Dependent = true )
        conditions      % (dependent) list of stim conditions
        
        tstat           % (dependent) t-stats of beta
        
        p               % (dependent) p-values of beta
        q               % (dependent) q-values of beta (false discovery rate)
        
    end
    
    methods
        % unique conditions
        function c = get.conditions( obj )
            if ~isempty(obj.variables)
                c = unique(obj.variables.cond); 
            else
                c = [];
            end
        end
        
        % t statistic calculation
        function tstat = get.tstat( obj )
            se=sqrt(sum(obj.covb_chol.^2,2));
            tstat = obj.beta ./ se;
            tstat(find(isnan(tstat)))=0;
        end
        
        % p value calculation
        function p = get.p( obj )
            t = obj.tstat;
            p = 2*tcdf(-abs(t), obj.dfe);
          
        end
        
        function pwr = power(obj,s)
            
            if(nargin<2)
                s='p<0.05';
            end
            
            tcrit=obj.getCritT(s);
            se=1/tcrit;
            t=1./sqrt(se.^2+obj.typeII_StdE.^2);
          
            pwr = 2*tcdf(-abs(t), obj.dfe);
        
        end
        
        % q values
        function q = get.q( obj )
            q = reshape( nirs.math.fdr( obj.p(:) )', size(obj.p) );
        end
        
        % critical value
        function out = getCritT( obj, s )
            %% getCritT - returns critical t-stat
            % 
            % Args:
            %     s - string specifying statistical significance
            %         (e.g. 'p < 0.05' or 'q < 0.1')
        
            s = strtrim( strsplit( s, '<' ) );
            
            if s{1} == 'p'
                pcrit = str2num(s{2});
                out = - tinv( abs(pcrit)/2, obj.dfe );
                
            elseif s{1} == 'q'
                t = abs(obj.tstat(:))';
                q = obj.q(:);
                
                [~,idx] = unique(q);
                
                out = interp1(q(idx), t(idx), str2num(s{2}));
            end
        end
        
        function out = table( obj )
            %% table - returns a table of the regression stats
            beta = obj.beta;
            tstat = obj.tstat;
            p = obj.p;
            q = obj.q;
            se = sqrt(sum(obj.covb_chol.^2,2));
            dfe = obj.dfe * ones(size(p));
            
            out = [obj.variables table(beta, se, tstat, dfe, p, q)];
        end
        
        function out = sorted( obj, colsToSortBy )
            %% sorted - returns sorted stats by columns in variables
            out = obj;
            if nargin < 2 | ~all(ismember(out.variables.Properties.VariableNames,colsToSortBy))
                colsToSortBy = {'VoxID', 'type', 'cond'};
            end
            
            [out.variables, idx] = sortrows(out.variables, colsToSortBy);
            
            out.beta = obj.beta(idx);
            out.covb_chol = obj.covb_chol(idx, :);
            out.typeII_StdE = obj.typeII_StdE(idx, :);
        end
        
        stats = ttest( obj, c, b );
        stats = ftest( obj, m );
        stats = jointTest( obj );
        
        h=draw( obj, vtype, vrange, thresh, powerthresh, viewpt );
        
        printAll( obj, vtype, vrange, thresh, powerthresh,viewpt, folder, ext );
    end
    
    methods (Access = protected)
        newNames = transformNames( obj, T );
    end
  
end