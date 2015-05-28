classdef MVGLM < nirs.modules.AbstractGLM
    properties
        useSpectralPriors = true;
        PPF = 50 / 5;
%         basis = Dictionary();
%         verbose
%         trend_func
    end
    
    methods
        function obj = MVGLM( prevJob )
            if nargin > 0, obj.prevJob = prevJob; end
            
            obj.name = 'Multivariate GLM';
            obj.basis('default') = nirs.design.basis.Canonical();
        end
        
        function S = runThis( obj, data )
            for i = 1:length(data)
       
                % get data
                d  = data(i).data;
                t  = data(i).time;
                Fs = data(i).Fs;
                
                % reshape data
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
                
                % distances
                l = data(i).probe.distances(idx);
                l = l(1:n:end);
                
                % new probe
                probe   = data(i).probe;
                lst     = link.type == link.type(1);
                probe.link = link(lst,:);
                probe.link.type = repmat( {'mv'}, [sum(lst) 1]);
                
                % outputs
                ncond = length(names);
                S(i) = nirs.ChannelStats();                   
                S(i).description    = data(i).description;
                S(i).names          = names';
                S(i).demographics   = data(i).demographics;
                S(i).probe          = probe;
                
                % fit data
                for iChan = 1:size(Y,3)
                    if obj.useSpectralPriors
                        thisX = X * l(iChan) * obj.PPF; 
                    else
                        thisX = X;
                    end
                    
                    stats = nirs.math.mv_ar_irls(thisX, Y(:,:,iChan), round(4*Fs), C);
                
                    % put stats
                    S(i).beta(:,iChan)      = stats.b(1:ncond);
                    S(i).covb(:,:,iChan)    = stats.covb(1:ncond, 1:ncond, :);
                    S(i).dfe                = stats.dfe;
                
                end
                
                % print progress
                obj.printProgress( i, length(data) )
            end

        end
    end
    
    methods ( Access = protected )
        function [X, names] = createX( obj, data )
            
            lambda = unique( data.probe.link.type );
            
            names = {};
            if obj.useSpectralPriors
                
                t       = data.time;
                stims   = data.stimulus;

                [xhbo, n] = nirs.design. ...
                    createDesignMatrix( stims, t, obj.basis, 'hbo' );
                
                names   = [names n];
                
                [xhbr, n] = nirs.design. ...
                    createDesignMatrix( stims, t, obj.basis, 'hbr' );
                
                names   = [names n];
                
                X = [];
                e = nirs.media.getspectra(lambda);
                
                for i = 1:length( lambda )
                    X(:,:,i) = [e(i,1)*xhbo e(i,2)*xhbr] * 1e-6;
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

