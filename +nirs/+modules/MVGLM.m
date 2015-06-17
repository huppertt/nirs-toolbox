classdef MVGLM < nirs.modules.AbstractGLM
    properties
        useSpectralPriors = true;
        PPF = 50 / 5;
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
                
                % sort data
                link = data(i).probe.link;
                [link, idx] = sortrows( link, {'source', 'detector', 'type'} );
                
                d = d(:,idx);
                
                % design mat
                [X, tbl] = nirs.design.createMultiVarDesignMat( ...
                    data(i).stimulus, t, obj.basis, link, true );
                
                C = obj.getTrendMatrix( t );
                
                % fit data
                %profile clear; profile on;
                stats = nirs.math.mv_ar_irls( X, d, round(4*Fs), C);
                %profile off; profile viewer
                
                
% %                 % unique sd pairs
% %                 [SD, ~, iSD] = unique(link(:,1:2), 'rows');
% %                 
% %                 % loop through sd pairs
% %                 X = []; names = {};
% %                 newLink = table([],[],[],[],'VariableNames',{'source', 'detector', 'type', 'cond'});
% %                 for j = 1:max(iSD)
% %                     lambda = link.type(iSD == j);
% %                     
% %                     [x, n, tbl] = obj.createX( data(i), lambda );
% %                     
% %                     X = blkdiag(X,x);
% %                     names = [names; n(:)];
% %                 end
% %                 
% %                 
% %                 
% %                 n = length( unique(link.type) );
% %                 
% %                 Y = [];
% %                 for j = 1:n
% %                    Y(:,j,:) = d(:,j:n:end); 
% %                 end
% %                 
% %                 % get experiment design
% %                 [X, names] = obj.createX( data(i) );
% %                 C = obj.getTrendMatrix( t );
% %                 
% %                 % distances
% %                 l = data(i).probe.distances(idx);
% %                 l = l(1:n:end);
                
                % new probe
                probe   = data(i).probe;
                lst     = link.type == link.type(1);
                probe.link = link(lst,:);
                probe.link.type = repmat( {'mv'}, [sum(lst) 1]);
                
                % outputs
                ncond = size(tbl,1);
                S(i) = nirs.core.ChannelStats();                   
                S(i).description    = data(i).description;
                S(i).variables          = tbl;
                S(i).demographics   = data(i).demographics;
                S(i).probe          = probe;
                S(i).beta = stats.b(1:ncond);
                S(i).covb = stats.covb(1:ncond, 1:ncond);
                S(i).dfe  = stats.dfe;
                
%                 % fit data
%                 for iChan = 1:size(Y,3)
%                     if obj.useSpectralPriors
%                         thisX = X * l(iChan) * obj.PPF; 
%                     else
%                         thisX = X;
%                     end
%                     
%                     stats = nirs.math.mv_ar_irls(thisX, Y(:,:,iChan), round(4*Fs), C);
%                 
%                     % put stats
%                     S(i).beta(:,iChan)      = stats.b(1:ncond);
%                     S(i).covb(:,:,iChan)    = stats.covb(1:ncond, 1:ncond, :);
%                     S(i).dfe                = stats.dfe;
%                 
%                 end
                
                % print progress
                obj.printProgress( i, length(data) )
            end

        end
    end
    
    methods ( Access = protected )
        function [X, names, link] = createX( obj, data, lambda )
          
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
        
        function createSpectralX( obj, data )
            
        end
        
        function createDodX( obj, data )
            
        end
    end
            
end

