classdef ChannelFStats
    %% CHANNELFSTATS - Holds F-statistics in channel space.
    % 
    % Properties: 
    %     description  - description of data (e.g. filename)
    %     variables    - a table containing source, detector, type, cond for
    %                    each F-stat
    %     F            - F statistics
    %     df1          - degrees of freedom 1
    %     df2          - degrees of freedom 2
    %     probe        - Probe object describing measurement geometry
    %     conditions   - (dependent) list of uqique conditions
    %     p            - (dependent) p-values of F
    %     q            - (dependent) q-values of F (false discovery rate)
    %     
    %  Methods:
    %     getCritF    - returns critical F value
    %     draw        - displays the time series data and stimulust timings
    %     table       - returns a table of all stats (minus full covariance)
    %     sorted      - returns new ChannelFStats object with sorted channels
    
    properties
        description
        variables	% description of data (e.g. filename)
        F           % F-stats
        df1         % degrees of freedom 1
        df2         % degrees of freedom 2
        probe       % probe geometry
         demographics    % Dictionary containing demographics info
    end
    
    properties ( Dependent = true )
        conditions  % (dependent) - a table containing source, detector, type, cond for each F-stat
        p           % (dependent) p-values of F
        q           % (dependent) q-values of F (false discovery rate)
    end
    
    methods
        % critical value
        function fcrit = getCritF( obj )
            %% getCritF - returns critical F-stat
            % 
            % Args:
            %     s - string specifying statistical significance
            %         (e.g. 'p < 0.05' or 'q < 0.1')
            fcrit = zeros(length(obj.names),1);
            for i = 1:length(obj.names)
                fcrit(i) = 1/finv( obj.pcrit, obj.df2(i), obj.df1(i) );
            end
        end
        
        % unique conditions
        function c = get.conditions( obj )
            if ~isempty(obj.variables)
                c = unique(obj.variables.cond); 
            else
                c = [];
            end
        end
        
        % p values
        function p = get.p( obj )
            p = zeros(size(obj.F));
            for i = 1:numel(obj.F)
                p(i) = fcdf( 1./obj.F(i), obj.df2(i), obj.df1(i) );
                
            end
        end
        
        % q values
        function q = get.q( obj )
            q = reshape( nirs.math.fdr( obj.p(:) )', size(obj.p) );
        end
        
        % stats table
        function out = table( obj )
            %% table - returns a table of all channelwise stats
            F = obj.F;
            p = obj.p;
            q = obj.q;
            df1 = obj.df1;
            df2 = obj.df2;
            
            out = [obj.variables table(F, df1, df2, p, q)];
        end
        
        % sorting
        function out = sorted( obj, colsToSortBy )
            %% sorted - returns sorted stats by columns in variables
            out = obj;
            if nargin < 2
                if(~ismember('source',out.variables.Properties.VariableNames) & ...
                        ismember('ROI',out.variables.Properties.VariableNames))
                    colsToSortBy = {'ROI', 'type', 'cond'};
                else
                    colsToSortBy = {'source', 'detector', 'type', 'cond'};
                end
            end
            
            [out.variables, idx] = nirs.util.sortrows(out.variables, colsToSortBy);
            
            out.F = obj.F(idx);
            out.df1 = obj.df1(idx);
            out.df2 = obj.df2(idx);
        end
        
        % draw
        f=draw( obj,  frange, thresh );
    end
end