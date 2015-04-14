classdef AbstractGLM < nirs.modules.AbstractModule
    %UNTITLED3 Summary of this class goes here
    %   Detailed explanation goes here
  
    properties
        basis       = HashTable();
        verbose     = false;
        isconstant  = true;
        trend_func  = @(t) nirs.design.trend.dctmtx(t, 1/125);
    end

    methods( Access = protected )
        
        % generate design matrix
        function [X, names] = createX( data )
            t       = data.time;
            stims   = data.stimulus;
            
            [X, names] = nirs.functional. ...
                createDesignMatrix( stims, t, obj.basis );
        end
        
        % generate baseline/trend regressors
        function C = getTrendMatrix( t )
            if ~isempty(obj.trend_func)
                C = obj.trend_func( t );
            else
                C = ones(size(t));
            end

            if obj.isconstant == false;
                C = C(:,2:end);
            end
        end
        
        % check rank and throw error
        function checkRank( X )
            if rank(X) < size(X,2)
                error( 'Design matrix is rank deficient.' )
            end
        end
                
        % check condition and issue warnings
        function checkCondtion( X )
            maxCond = 100; 
            if cond([X C]) > maxCond
                warning(['high collinearity: cond(X) = ' num2str(cond([X C]))])
            end
        end
        
        % print progress
        function printProgress(n, N)
            if obj.verbose
                fprintf( 'Finished %4i of %4i.\n', n, N )
            end
        end        
    end
    
end

