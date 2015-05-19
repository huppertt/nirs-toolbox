classdef MVGLM < nirs.modules.AbstractGLM
    properties
        useSpectralPriors = true;
    end
    
    methods
        function obj = MVGLM( prevJob )
            if nargin > 0, obj.prevJob = prevJob; end
            
            obj.name = 'Multivariate GLM';
            obj.basis('default') = nirs.design.basis.Canonical();
        end
        
        function S = runThis( obj, data )
            for i = 1:length(data)
                
                %% GET DATA AND ORGANIZE INTO 3D ARRAY
                %% loop through channels and call mv_ar_irls
                %% override createX
                
                % get data
                d  = data(i).data;
                t  = data(i).time;
                Fs = data(i).Fs;
                
                % get experiment design
                [X, names] = obj.createX( data(i) );
                C = obj.getTrendMatrix( t );
                
                % check model
                obj.checkRank( [X C] )
                obj.checkCondition( [X C] )
                
                % run regression
                thisS = nirs.math.ar_irls( d, [X C], round(4*Fs) );
                
                % put stats
                ncond = length(names);
                S(i) = nirs.ChannelStats();
                S(i).description = data(i).description;
                S(i).beta = thisS.beta(1:ncond,:);
                S(i).covb = thisS.covb(1:ncond, 1:ncond, :);
                S(i).dfe  = thisS.dfe(1);
                
                S(i).names          = names';
                S(i).demographics   = data(i).demographics;
                S(i).probe          = data(i).probe;
                
                % print progress
                obj.printProgress( i, length(data) )
            end

        end
        
    end
    
end

