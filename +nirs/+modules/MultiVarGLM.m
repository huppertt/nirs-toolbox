classdef MultiVarGLM < nirs.modules.AbstractGLM

    methods
        function obj = MultiVarGLM( prevJob )
            if nargin > 0, obj.prevJob = prevJob; end

            obj.name             	= 'Multivariate GLM';
            obj.basis('default')    = nirs.design.basis.Canonical();
            obj.trend_func          = @nirs.design.trend.constant;
        end

        function S = runThis( obj, data )
            vec = @(x) x(:);

            for i = 1:numel(data)

                % get data
                d  = data(i).data;
                t  = data(i).time;
                Fs = data(i).Fs;

                % sort data
                link = data(i).probe.link;
                if(~isempty(strfind(class(data(i).probe),'nirs')))
                    if(~ismember('source',link.Properties.VariableNames) & ...
                            ismember('ROI',link.Properties.VariableNames))
                        [link, idx] = nirs.util.sortrows(link, {'ROI','type'});
                    else
                        [link, idx] = nirs.util.sortrows(link, {'source', 'detector','type'});
                    end
                elseif(~isempty(strfind(class(probe),'eeg')))
                    [link, idx] = nirs.util.sortrows(link, {'electrode','type'});
                else
                    error('data type not supported');
                end

                d = d(:,idx);

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
                    % fit data
                    stats = nirs.math.mv_ar_irls( U(:,lst)*s(lst,lst), d, round(4*Fs));
                    stats.beta=V*stats.beta;
                    IV=kron(V,eye(size(d,2)));
                    stats.covb=IV*stats.covb*IV';
                else

                    % fit data
                    stats = nirs.math.mv_ar_irls( [X C], d, round(4*Fs));
                end

                if(~isempty(V1))
                    n=size(V1,2);
                    stats.beta=V1*stats.beta(1:n,:);
                    IV=kron(V1,eye(size(d,2)));
                    stats.covb=IV*stats.covb*IV';
                end

                probe=data(i).probe;
                % put stats
                ncond = length(names);
                nchan = size(probe.link, 1);

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

                %                 if(obj.useFstats)
                %                     S(i).pvalue_fixed=vec(stats.Fpval(1:ncond,:)' );
                %                 end

                
%                 lst=[];
%                 for i=1:ncond
%                     lst=[lst i:size(stats.beta,1):size(stats.covb,1)];
%                 end
%                 lst=sort(lst);
                lst=[1:length(S(i).beta)];

                covb = stats.covb(lst,lst);
                
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

    end

end