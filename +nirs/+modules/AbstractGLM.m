classdef AbstractGLM < nirs.modules.AbstractModule
%% ABSTRACTGLM  -  Abstract class with common GLM functions.
  
    properties
        basis       = Dictionary();                 % dictionary of bases using stim name as key
        verbose     = true;                         % flag to show GLM progress
        trend_func  = @nirs.design.trend.constant   % a function of time that returns trend regressors
        goforit = false;  % flag to force processing of poorly conditioned design matrices (use with caution)
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
        end
        
        % print progress
        function printProgress(obj, n, N)
            if obj.verbose
                nirs.util.flushstdout(1);
                fprintf( 'Finished %4i of %4i.\n', n, N )
            end
        end  
    end
    
    methods
        % check rank and throw error
        function checkRank(obj,X )
            if rank(X) < size(X,2)
                if(~obj.goforit)
                    error( 'Design matrix is rank deficient.' )
                else
                    disp( 'Design matrix is rank deficient.' )
                end
            end
        end
                
        % check condition and issue warnings
        function checkCondition(obj,X )
            maxCond = 300; 
            if cond(X) > maxCond
                disp(['High collinearity: cond(X) = ' num2str(cond(X)) '.'])
            end
        end
    end
    
end

