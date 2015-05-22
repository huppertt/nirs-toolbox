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
                
                link = data(i).probe.link;
                [link, idx] = sortrows( link, {'source', 'detector'} );
                
                d = d(:,idx);
                
                n = length( unique(link.type) );
                
                Y = [];
                for j = 1:n
                   Y(:,j,:) = d(:,j:n:end); 
                end
                
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
    
    methods ( Access = protected )
        function [X, names] = createX( obj, data )
            
            X = []; names = {};
            if obj.useSpectralPriors
                
                types = {'hbo','hbr'};
                for i = 1:length(types)
                    t       = data.time;
                    stims   = data.stimulus;

                    [x, n] = nirs.design. ...
                        createDesignMatrix( stims, t, obj.basis, types{i} );

                    X       = [X x];
                    names   = [names n];
                end
            else
                t       = data.time;
                stims   = data.stimulus;
                
                [X, names] = nirs.design. ...
                        createDesignMatrix( stims, t, obj.basis );
            end
            
        end
    end
            
end

