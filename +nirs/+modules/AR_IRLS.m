classdef AR_IRLS < nirs.modules.AbstractGLM
   
    methods
        function obj = AR_IRLS( prevJob )
            if nargin > 0, obj.prevJob = prevJob; end
            
            obj.name = 'GLM via AR(P)-IRLS';
            obj.basis('default') = nirs.design.basis.Canonical();
        end
        
        function S = runThis( obj, data )
            vec = @(x) x(:);
            
            for i = 1:length(data)
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
                stats = nirs.math.ar_irls( d, [X C], round(4*Fs) );
                
                % put stats
                ncond = length(names);
                nchan = size(data(i).probe.link, 1);
                
                link = repmat( data(i).probe.link, [ncond 1] );
                condition = repmat(names, [nchan 1]);
                                
                S(i) = nirs.core.ChannelStats();
                
                S(i).variables = [link table(condition)];
                S(i).beta = vec( stats.beta(1:ncond,:) );
                
                covb = [];
                for j = 1:nchan
                   covb = blkdiag(covb, stats.covb(1:ncond, 1:ncond, j)); 
                end
                
                S(i).covb = covb;
                
%                 S(i).covb = thisS.covb(1:ncond, 1:ncond, :);
                S(i).dfe  = stats.dfe(1);
                
                S(i).description = data(i).description;
                
%                 S(i).names          = names';
                S(i).demographics   = data(i).demographics;
                S(i).probe          = data(i).probe;
                
                % print progress
                obj.printProgress( i, length(data) )
            end

        end
        
    end
    
end

