
classdef AR_IRLS < nirs.modules.AbstractGLM
    %% AR_IRLS - Performs first-level per file GLM analysis.
    %
    % Options:
    %     basis       - a Dictionary object containing temporal bases using stim name as key
    %     verbose     - flag to display progress
    %     trend_func  - a function that takes in a time vector and returns trend regressors
    %
    % Example:
    %     j = nirs.modules.AR_IRLS();
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
    properties(Hidden = true)
        useREML=false;
        nonstationary_noise=false;
        useFstats=false;
        useGPU=false;
        precisionSingle=false;
        Pmax='4*Fs';
    end
    methods
        function obj = AR_IRLS( prevJob )
            if nargin > 0, obj.prevJob = prevJob; end
            
            obj.name = 'GLM via AR(P)-IRLS';
            obj.basis('default') = nirs.design.basis.Canonical();
            obj.citation=['Barker, Jeffrey W., Ardalan Aarabi, and Theodore J. Huppert.'...
                '"Autoregressive model based algorithm for correcting motion and serially '...
                'correlated errors in fNIRS." Biomedical optics express 4.8 (2013): 1366-1379.'];
            
        end
        
        function S = runThis( obj, data )
            vec = @(x) x(:);
            
            for i = 1:numel(data)
                % get data
                d  = data(i).data;
                d=d-ones(size(d,1),1)*nanmean(d,1);
                d=d-ones(size(d,1),1)*nanmean(d,1);
                t  = data(i).time;
                Fs = data(i).Fs;
                
                probe = data(i).probe;
                
                Pmax = obj.Pmax;
                if(isstr(Pmax))
                    eval(['Pmax=' Pmax ';']);
                end
                Pmax=ceil(Pmax);    
                
                % make sure data is in order
                
                if(~isempty(strfind(class(probe),'nirs')))
                    if(~ismember('source',probe.link.Properties.VariableNames) & ...
                            ismember('ROI',probe.link.Properties.VariableNames))
                        [probe.link, idx] = nirs.util.sortrows(probe.link, {'ROI','type'});
                    else
                        [probe.link, idx] = nirs.util.sortrows(probe.link, {'source', 'detector','type'});
                    end
                elseif(~isempty(strfind(class(probe),'eeg')))
                    [probe.link, idx] = nirs.util.sortrows(probe.link, {'electrode','type'});
                else
                    error('data type not supported');
                end
                d = d(:, idx);
                
                % get experiment design
                [X, names] = obj.createX( data(i) );
                C = obj.getTrendMatrix( t );
                
                X(find(isnan(X)))=0;
                lstNull=find(all(X==0,1));
                if(~isempty(lstNull))
                    for id=1:length(lstNull)
                        disp(['Stim Condition : ' names{lstNull(id)} ' is all zeros- removing']);
                    end
                    X(:,lstNull)=[];
                    lstA=1:length(names); lstA(lstNull)=[];
                    names={names{lstA}};
                end
                
                if(cond(X)>300 & obj.goforit)
                    [U1,S1,V1]=nirs.math.mysvd(X);
                    lst=1:rank(X);
                    X = U1(:,lst)*S1(lst,lst);
                    V1=V1(:,lst);
                else
                    V1=[];
                    % check model
                    obj.checkRank( [X C] )
                    obj.checkCondition( [X C] )
                end
                
                
                
                if(rank([X C]) < size([X C],2) & obj.goforit)
                    disp('Using PCA regression model');
                    [U,s,V]=nirs.math.mysvd([X C]);
                    lst=find(diag(s)>eps(1)*100);
                    V=V(:,lst);
                    if(obj.useREML)
                        stats = nirs.math.ar_irls_REML( d, U(:,lst)*s(lst,lst), Pmax ,[],obj.useGPU,obj.precisionSingle);
                    else
                        if(obj.nonstationary_noise)
                            stats = nirs.math.ar_irnnls( d, U(:,lst)*s(lst,lst), Pmax ,[], obj.useGPU, obj.precisionSingle );
                        else
                            stats = nirs.math.ar_irls( d, U(:,lst)*s(lst,lst), Pmax ,[],false,obj.useGPU, obj.precisionSingle);
                        end
                    end
                    stats.beta=V*stats.beta;
                    c=[];
                    for j=1:size(stats.covb,3)
                        for j2=1:size(stats.covb,4)
                            c(:,:,j,j2)=V*squeeze(stats.covb(:,:,j,j2))*V';
                        end
                    end
                    stats.covb=c;
                else
                    
                    % run regression
                    if(obj.useREML)
                        stats = nirs.math.ar_irls_REML( d, [X C], Pmax,[],obj.useGPU, obj.precisionSingle );
                    else
                        if(obj.nonstationary_noise)
                             stats = nirs.math.ar_irnsls( d, [X C], Pmax ,[],obj.useGPU, obj.precisionSingle);
                        else
                            if(obj.useFstats)
                                stats = nirs.math.ar_irls_ftest( d, [X C], Pmax ,[],obj.useGPU, obj.precisionSingle);
                            else
                                stats = nirs.math.ar_irls( d, [X C], Pmax ,[],false,obj.useGPU, obj.precisionSingle);
                            end
                        end
                    end
                end
                
                
                
                if(~isempty(V1))
                    n=size(V1,2);
                    stats.beta=V1*stats.beta(1:n,:);
                    c=[];
                    for j=1:size(stats.covb,3)
                        for j2=1:size(stats.covb,4)
                            c(:,:,j,j2)=V1*squeeze(stats.covb(1:n,1:n,j,j2))*V1';
                        end
                    end
                    stats.covb=c;
                end
                
                % put stats
                ncond = length(names);
                nchan = size(data(i).probe.link, 1);
                
                link = repmat( probe.link, [ncond 1] );
                condition = repmat(names(:)', [nchan 1]);
                condition = condition(:);
                
                if(~isempty(strfind(class(probe),'nirs')))
                    S(i) = nirs.core.ChannelStats();
                elseif(~isempty(strfind(class(probe),'eeg')))
                    S(i) = eeg.core.ChannelStats();
                else
                    warning('unsupported data type');
                    S(i) = nirs.core.ChannelStats();
                end
                S(i).variables = [link table(condition,'VariableNames',{'cond'})];
                S(i).beta = vec( stats.beta(1:ncond,:)' );
                
                if(obj.useFstats)
                    S(i).pvalue_fixed=vec(stats.Fpval(1:ncond,:)' );
                end
                
                covb = zeros( nchan*ncond );
                if(ndims(stats.covb)==4)
                    for j = 1:nchan
                        for k=1:nchan
                            idx = (0:ncond-1)*nchan + j;
                            idx2 = (0:ncond-1)*nchan + k;
                            covb(idx, idx2) = stats.covb(1:ncond, 1:ncond, j,k);
                        end
                    end
                elseif(ndims(stats.covb)==3)
                    for j = 1:nchan
                        idx = (0:ncond-1)*nchan + j;
                        covb(idx, idx) = stats.covb(1:ncond, 1:ncond, j);
                        
                    end
                else
                    covb = stats.covb(1:ncond*nchan,1:ncond*nchan);
                end
                
                %ensure positive/definant (sometimes off due to numerical
                %prec.
                
                lst=find(~all(isnan(covb),1));
                [U,s,V]=svd(covb(lst,lst));
                covb(lst,lst)=0.5*(U*s*U'+V*s*V');
                
                S(i).covb = covb;
                
                S(i).dfe  = stats.dfe(1);
                
                S(i).description = data(i).description;
                
                S(i).demographics   = data(i).demographics;
                S(i).probe          = probe;
                
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
                if(obj.verbose)
                    obj.printProgress( i, length(data) )
                end
            end
            
        end
        
        
        function prop = javaoptions(obj)
            
            prop=javaoptions@nirs.modules.AbstractGLM(obj);
            opts=obj.options;
            
            diction=nirs.util.createDictionaryFromToolset('nirs.design.basis');
            DictionaryProp=javatypes('enum',{diction.values});
            set(DictionaryProp,'Name','basis','Value','test');
            set(DictionaryProp,'Category','Misc');
            set(DictionaryProp,'Description','Select the canonical basic function');
            prop(find(ismember(opts,'basis')))=DictionaryProp;
            
            
            
            
        end
        
    end
    
end

