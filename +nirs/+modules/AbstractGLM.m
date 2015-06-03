classdef AbstractGLM < nirs.modules.AbstractModule
    %UNTITLED3 Summary of this class goes here
    %   Detailed explanation goes here
  
    properties
        basis       = Dictionary();
        verbose     = false;
        isconstant  = true;
        trend_func  = @(t) nirs.design.trend.legendre(t, 3);
    end

    methods( Access = protected )
        
        % generate design matrix
        function [X, names] = createX( obj, data )
            t       = data.time;
            stims   = data.stimulus;
            
            [X, names] = nirs.design. ...
                createDesignMatrix( stims, t, obj.basis );
        end
        
        % generate baseline/trend regressors
        function C = getTrendMatrix( obj, t )
            if ~isempty(obj.trend_func)
                C = obj.trend_func( t );
            else
                C = ones(size(t));
            end

            if obj.isconstant == false;
                C = C(:,2:end);
            end
        end
        
        % print progress
        function printProgress(obj, n, N)
            if obj.verbose
                fprintf( 'Finished %4i of %4i.\n', n, N )
            end
        end  
    end
    
    methods( Access = protected, Static )
        % check rank and throw error
        function checkRank( X )
            if rank(X) < size(X,2)
                error( 'Design matrix is rank deficient.' )
            end
        end
                
        % check condition and issue warnings
        function checkCondition( X )
            maxCond = 100; 
            if cond(X) > maxCond
                warning(['High collinearity: cond(X) = ' num2str(cond(X)) '.'])
            end
        end
    end
    
end

