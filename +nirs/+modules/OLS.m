classdef OLS < nirs.modules.AbstractGLM
%% OLS - Performs first-level per file GLM analysis via Ordinary Least Squares.
%
% This is surely not an appropriate GLM analysis and is only here for ROC
% testing. 
%
% Options:
%     basis       - a Dictionary object containing temporal bases using stim name as key
%     verbose     - flag to display progress
%     trend_func  - a function that takes in a time vector and returns trend regressors
%     
% Example:
%     j = nirs.modules.OLS();
%     
%     b = Dictionary();
%     b('default') = nirs.design.basis.Canonical(); % default basis
%     b('A')       = nirs.design.basis.Gamma();     % a different basis for condition 'A'
%     
%     j.basis = b;
%     
%     j.trend_func = @(t) nirs.design.trend.legendre(t, 3); % 3rd order polynomial for trend
%     
% Note: 
%     trend_func must at least return a constant term unless all baseline periods are
%     specified explicitly in the stimulus design with BoxCar basis functions
    
    methods
        function obj = OLS( prevJob )
            if nargin > 0, obj.prevJob = prevJob; end
            
            obj.name = 'GLM via OLS';
            obj.basis('default') = nirs.design.basis.Canonical();
        end
        
        function S = runThis( obj, data )
            vec = @(x) x(:);
            
            for i = 1:length(data)
                % get data
                d  = data(i).data;
                t  = data(i).time;
                Fs = data(i).Fs;
                
                probe = data(i).probe;
                
                % make sure data is in order
                [probe.link, idx] = sortrows(probe.link, {'source', 'detector', 'type'});
                d = d(:, idx);
                
                % get experiment design
                [X, names] = obj.createX( data(i) );
                C = obj.getTrendMatrix( t );
                
                % check model
                obj.checkRank( [X C] )
                obj.checkCondition( [X C] )
                
                % run regression
                %stats = nirs.math.ar_irls( d, [X C], round(4*Fs) );
                stats.beta = [X C]\ d;
                for j = 1:size(d,2)
                    stats.covb(:,:,j) = pinv([X C]'*[X C]) * var(d(:,j) - [X C]*stats.beta(:,j));
                end
                stats.dfe = size(d,1) - rank([X C]);
                
                
                % put stats
                ncond = length(names);
                nchan = size(data(i).probe.link, 1);
                
                link = repmat( probe.link, [ncond 1] );
                cond = repmat(names(:)', [nchan 1]);
                cond = cond(:);
                                                
                S(i) = nirs.core.ChannelStats();
                
                S(i).variables = [link table(cond)];
                S(i).beta = vec( stats.beta(1:ncond,:)' );
                
                covb = zeros( nchan*ncond );
                for j = 1:nchan
                   idx = (0:ncond-1)*nchan + j;
                   covb(idx, idx) = stats.covb(1:ncond, 1:ncond, j);
                end
                
                S(i).covb = covb;
                
                S(i).dfe  = stats.dfe(1);
                
                S(i).description = data(i).description;
                
                S(i).demographics   = data(i).demographics;
                S(i).probe          = data(i).probe;
                
                % print progress
                obj.printProgress( i, length(data) )
            end

        end
        
    end
    
end

