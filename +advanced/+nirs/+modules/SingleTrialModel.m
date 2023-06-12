classdef SingleTrialModel < nirs.modules.AbstractGLM
%% Performs single trial first-level per file GLM analysis using a random effects model.
% 
% Options:
%     basis       - a Dictionary object containing temporal bases using stim name as key
%     verbose     - flag to display progress
%     trend_func  - a function that takes in a time vector and returns trend regressors
%     
% Example:
%     j = nirs.modules.SingleTrialModel();
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
    properties
        PCA_reduction = 0;
        randomeffectsmodel=false;
    end
    methods
        function obj = SingleTrialModel( prevJob )
            if nargin > 0, obj.prevJob = prevJob; end
            
            obj.name = 'Single Trial GLM Model';
            obj.basis('default') = nirs.design.basis.Canonical();
            obj.goforit=true;
            
        end
        
        function S = runThis( obj, data )
            vec = @(x) x(:);
            
            for i = 1:numel(data)
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
                    error('unsupported data type');
                end
                d = d(:, idx);
                
                 % get experiment design
                [X, names] = obj.createX( data(i) );
                C = obj.getTrendMatrix( t );
                
                X(find(isnan(X)))=0;
                
                % split the conditions to create the random effects model
                keys=data(i).stimulus.keys;
                if(obj.randomeffectsmodel)
                        stimulus=data(i).stimulus;
                    else
                        stimulus=Dictionary;
                end

                ttests={}; cnt=1;
                for k=1:length(keys)
                    name=keys{k};
                    
                    %if including multiple terms per condition
                    str2={};
                    for id=1:length(names)
                        if(~isempty(strfind(names{id},[name ':'])))
                            str2{end+1}=names{id}(strfind(names{id},[name ':'])+length([name]):end);
                        end
                    end
                    if(isempty(str2))
                        str2{1}='';
                    end
                    
                    stim=data(i).stimulus(name);
                    disp(['Split ' name ' into ' num2str(length(stim.onset)) ' trials']);
                    
                    
                    for l=1:length(stim.onset)
                        stim=data(i).stimulus(name);
                        stim.dur([1:l-1 l+1:end])=[];
                        stim.amp([1:l-1 l+1:end])=[];
                        stim.onset([1:l-1 l+1:end])=[];
                        stimulus([name '_trial' num2str(l)])=stim;
                        for id=1:length(str2)
                            if(obj.randomeffectsmodel)
                                ttests{cnt,1}=[name str2{id} '+' name '_trial' num2str(l) str2{id}];
                                ttests{cnt,2}=[name '_trial' num2str(l) str2{id}];
                            else
                                ttests{cnt,1}=[name '_trial' num2str(l) str2{id}];
                                ttests{cnt,2}=[name '_trial' num2str(l) str2{id}];
                            end
                            
                            cnt=cnt+1;
                        end
                    end
                    
                end
                data(i).stimulus=stimulus;
                
                % get experiment design
                [X, names] = obj.createX( data(i) );
                C = obj.getTrendMatrix( t );
                
                X(find(isnan(X)))=0;
               
                
                % run regression
                if(obj.PCA_reduction>0 || (rank([X C]) < size([X C],2) & obj.goforit))
                    disp('Using PCA regression model');
                    [U,s,V]=nirs.math.mysvd([X C]);
                    lst=nirs.math.sig_eigens(s,0.05);
                    lst((1-cumsum(diag(s(1:length(lst),1:length(lst)))./sum(diag(s))))<obj.PCA_reduction)=[];
                    V=V(:,lst);
                    
                    stats = nirs.math.ar_irls( d, U(:,lst)*s(lst,lst), round(4*Fs) );
                    stats.beta=V*stats.beta;
                    for j=1:size(stats.covb,3)
                        c(:,:,j)=V*squeeze(stats.covb(:,:,j))*V';
                    end
                    stats.covb=c;
                else
                    %check model
                    obj.checkRank( [X C] )
                    obj.checkCondition( [X C] )
                    % run regression
                    stats = nirs.math.ar_irls( d, [X C], round(4*Fs) );
                end
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
                
                [U,s,V]=svd(covb);
                covb=0.5*(U*s*U'+V*s*V');
                
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
                
                
                % Put the model back together
                S(i)=S(i).ttest({ttests{:,1}},[],{ttests{:,2}});
                % print progress
                obj.printProgress( i, length(data) )
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

