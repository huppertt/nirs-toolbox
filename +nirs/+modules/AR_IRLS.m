classdef AR_IRLS < nirs.modules.AbstractGLM
   
    methods
        function obj = AR_IRLS( prevJob )
           obj.name = 'GLM via AR(P)-IRLS';
           obj.basis('default') = nirs.design.basis.Canonical();
           
           if nargin > 0
               obj.prevJob = prevJob;
           end
        end
        
        function stats = runThis( obj, data )
            for i = 1:length(data)
                
                d  = data(i).data;
                t  = data(i).time;
                Fs = data(i).Fs;
                
                [X, names] = createX( data );
                
                C = getTrendMatrix( t );
                
                checkRank( [X C] )
                
                checkCondition( [X C] )
                
                warning('off','stats:statrobustfit:IterationLimit')
                thisS = nirs.math.ar_irls( d, [X C], round(4*Fs) );
                thisS.X = X;
                thisS.C = C;
                thisS.names = names';
                thisS.stimulus = data(i).stimulus;
                thisS.demographics = data(i).demographics;
                thisS.probe = data(i).probe;
                
                S(i) = thisS;
                
                printProgress( i, length(data) )
            end
            
            % output stats
            stats = S;
                
%                 % generate design matrix
%                 stims = data(i).stimulus;
%                 [X, names] = nirs.functional. ...
%                     createDesignMatrix( stims, t, obj.basis );
%                 
%                 % generate baseline/trend regressors
%                 if ~isempty(obj.trend_func)
%                     C = obj.trend_func( t );
%                 else
%                     C = ones(size(t));
%                 end
%                 
%                 if obj.isconstant == false;
%                     C = C(:,2:end);
%                 end
%                 
%                 % check rank
%                 if rank(X) < size(X,2)
%                     error( 'Design matrix is rank deficient.' )
%                 end
%                 
%                 % check condition
%                 maxCond = 100; 
%                 if cond([X C]) > maxCond
%                     warning(['high collinearity: cond(X) = ' num2str(cond([X C]))])
%                 end
%                                 
%                 % call ar_irls
%                 warning('off','stats:statrobustfit:IterationLimit')
%                 thisS = ar_irls( d, [X C], round(4*Fs) );
%                 thisS.X = X;
%                 thisS.C = C;
%                 thisS.names = names';
%                 thisS.stimulus = data(i).stimulus;
%                 thisS.demographics = data(i).demographics;
%                 thisS.probe = data(i).probe;
%                 
%                 S(i) = thisS;
%                 
%                 % output stats
%                 stats = S;
%                 
%                 if obj.verbose
%                     fprintf( 'Finished %4i of %4i.\n', i, length(data) )
%                 end
% 
%             end
        end
        
    end
    
end

