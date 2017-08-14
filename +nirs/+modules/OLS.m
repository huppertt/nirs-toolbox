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
            obj.citation='Huppert, T. J., Diamond, S. G., Franceschini, M. A., & Boas, D. A. (2009). HomER: a review of time-series analysis methods for near-infrared spectroscopy of the brain. Applied optics, 48(10), D280-D298.';
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
                if(~isempty(strfind(class(probe),'nirs')))
                    [probe.link, idx] = sortrows(probe.link, {'source', 'detector','type'});
                elseif(~isempty(strfind(class(probe),'eeg')))
                    [probe.link, idx] = sortrows(probe.link, {'electrode','type'});
                else
                    error('data type not supported');
                end
                d = d(:, idx);
                
                % get experiment design
                [X, names] = obj.createX( data(i) );
                C = obj.getTrendMatrix( t );
                
                % check model
                obj.checkRank( [X C] )
                obj.checkCondition( [X C] )
                
                % run regression
                if(rank([X C]) < size([X C],2) & obj.goforit)
                    disp('Using PCA regression model');
                    [U,s,V]=nirs.math.mysvd([X C]);
                    lst=find(diag(s)>eps(1)*10);
                    V=V(:,lst);
                    A=U(:,lst)*s(lst,lst);
                    stats.beta = A\ d;
                    for j = 1:size(d,2)
                        stats.covb(:,:,j) = pinv(A'*A) * var(d(:,j) - A*stats.beta(:,j));
                    end
                    
                    stats.beta=V*stats.beta;
                    for j=1:size(stats.covb,3)
                        c(:,:,j)=V*squeeze(stats.covb(:,:,j))*V';
                    end
                    stats.covb=c;
                    stats.dfe = size(d,1) - rank(A);
                else
                    stats.beta = [X C]\ d;
                    for j = 1:size(d,2)
                        stats.covb(:,:,j) = pinv([X C]'*[X C]) * var(d(:,j) - [X C]*stats.beta(:,j));
                    end
                    stats.dfe = size(d,1) - rank([X C]);
                end
                
                %stats = nirs.math.ar_irls( d, [X C], round(4*Fs) );

                
                
                
                % put stats
                ncond = length(names);
                nchan = size(data(i).probe.link, 1);
                
                link = repmat( probe.link, [ncond 1] );
                cond = repmat(names(:)', [nchan 1]);
                cond = cond(:);
                                                
                if(~isempty(strfind(class(probe),'nirs')))
                    S(i) = nirs.core.ChannelStats();
                elseif(~isempty(strfind(class(probe),'eeg')))
                    S(i) = eeg.core.ChannelStats();
                else
                    warning('unsupported data type');
                    S(i) = nirs.core.ChannelStats();
                end
                
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
                
                stim=Dictionary;
                for j=1:data(i).stimulus.count;
                    ss=data(i).stimulus.values{j};
                    if(isa(ss,'nirs.design.StimulusEvents'))
                        s=nirs.design.StimulusEvents;
                        s.name=ss.name;
                        s.dur=mean(ss.dur);
                        stim(data(i).stimulus.keys{j})=s;
                    end
                end
                
                S(i).basis.base=obj.basis;
                S(i).basis.Fs=Fs;
                S(i).basis.stim=stim;
                
                
                
                % print progress
                obj.printProgress( i, length(data) )
            end

        end
        
    end
    
end

