classdef AR_IRLS < nirs.functional.AbstractModule
    %UNTITLED3 Summary of this class goes here
    %   Detailed explanation goes here
  
    properties
        
        basis       = nirs.HashTable;
        constant    = true;
        verbose     = false;
        % hpf_Fc      = 1/125;
        trend_func  = @(t) nirs.functional.dctmtx(t, 1/125);
    end
    
    methods

        function obj = AR_IRLS( prevJob )
           obj.name = 'GLM via AR(P)-IRLS';
           obj.basis('default') = nirs.functional.basis.Canonical();
           
           if nargin > 0
               obj.prevJob = prevJob;
           end
        end
        
        function stats = execute( obj, data )
            for i = 1:length(data)
                
                d = data(i).data;
                t = data(i).time;
                Fs = data(i).Fs;
                
                % generate design matrix
                stims = data(i).stimulus;
                [X, names] = nirs.functional. ...
                    createDesignMatrix( stims, t, obj.basis );
                
                % generate baseline/trend regressors
%                 if strcmpi('trend','dct')
%                     C = nirs.functional.dctmtx( t, obj.trend_param )';
%                 elseif stcmpi('trend','poly')
%                     C = nirs.functional.legendre( t, obj.trend_param )';
%                 end
                C = obj.trend_func( t );
                
                if obj.constant == false;
                    C = C(:,2:end);
                end
                
%                 for j = 1:size(C,2)
%                     names{end+1} = ['dct_' num2str(j)];
%                 end
                
                % check rank
                if rank(X) < size(X,2)
                    error( 'Design matrix is rank deficient.' )
                end
                
                % check condition
                maxCond = 30; 
                if cond([X C]) > maxCond
                    warning('Consider lowering HPF cutoff to reduce collinearity.')
                end
                                
                % call ar_irls
                warning('off','stats:statrobustfit:IterationLimit')
                thisS = ar_irls( d, [X C], round(4*Fs) );
                thisS.X = X;
                thisS.C = C;
                thisS.names = names';
                thisS.stimulus = data(i).stimulus;
                thisS.demographics = data(i).demographics;
                thisS.probe = data(i).probe;
                
                S(i) = thisS;
                
                % output stats
                stats = S;
                
                if obj.verbose
                    fprintf( 'Finished %4i of %4i.\n', i, length(data) )
                end
            end
        end
        
        function options = getOptions( obj )
            options = [];
        end
           
        function obj = putOptions( obj, options )
        end
        
    end
    
end

