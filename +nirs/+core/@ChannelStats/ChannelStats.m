classdef ChannelStats
    %% CHANNELSTATS - Holds regression stats for channel space.
    % 
    % Properties: 
    %     description  - description of data (e.g. filename)
    %     variables    - a table containing w/ source, detector, type, cond for
    %                    each regression coefficient
    %     beta         - regression coefficients
    %     covb         - covariance of beta
    %     dfe          - degrees of freedom
    %     probe        - Probe object describing measurement geometry
    %     demographics - Dictionary containing demographics info
    %                    (e.g. demographics('age') returns 28)
    %     conditions   - (dependent) list of stim conditions
    %     tstat        - (dependent) t-stats of beta
    %     p            - (dependent) p-values of beta
    %     q            - (dependent) q-values of beta (false discovery rate)
    %     
    %  Methods:
    %     getCritT    - returns critical t value
    %     draw        - displays the time series data and stimulust timings
    %     table       - returns a table of all stats (minus full covariance)
    %     ttest       - performs a t-test and returns a new ChannelStats object
    %     ftest       - performs an F-test across conditions 
    %     jointTest   - performs a joint hypothesis test across all channels in
    %                   each SD pair (i.e. joint test of hbo & hbr for S1-D1)

    properties
        description     % description of data (e.g. filename)
        
        variables       % table describing regression coefficients
        
        beta            % regression coefficients
        covb            % covariance of beta
        dfe          	% degrees of freedom
        
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
            tstat = obj.beta ./ sqrt(diag(obj.covb));
        end
        
        % p value calculation
        function p = get.p( obj )
            t = obj.tstat;
            p = 2*tcdf(-abs(t), obj.dfe);
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
            se = sqrt(diag(obj.covb));
            dfe = obj.dfe * ones(size(p));
            
            out = [obj.variables table(beta, se, tstat, dfe, p, q)];
        end
        
        stats = ttest( obj, c, b );
        stats = ftest( obj, m );
        stats = jointTest( obj );
        
        draw( obj, vtype, vrange, thresh );
    end
    
    methods (Access = protected)
        newNames = transformNames( obj, T );
    end
  
end