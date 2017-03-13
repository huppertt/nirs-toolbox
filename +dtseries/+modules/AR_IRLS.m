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
                
                W = chol(data(i).cov);
                d=d*inv(W)';
                
                % get experiment design
                [X, names] = obj.createX( data(i) );
                C = obj.getTrendMatrix( t );
                
                X(find(isnan(X)))=0;
                
                % check model
                obj.checkRank( [X C] )
                obj.checkCondition( [X C] )
                
                if(rank([X C]) < size([X C],2) & obj.goforit)
                    disp('Using PCA regression model');
                    [U,s,V]=nirs.math.mysvd([X C]);
                    lst=find(diag(s)>eps(1)*10);
                    V=V(:,lst);
                    stats = nirs.math.ar_irls( d, U(:,lst)*s(lst,lst), round(4*Fs) );
                    stats.beta=V*stats.beta;
                    for j=1:size(stats.covb,3)
                        c(:,:,j)=V*squeeze(stats.covb(:,:,j))*V';
                    end
                    stats.covb=c;
                else
                    
                    % run regression
                    stats = nirs.math.ar_irls( d, [X C], round(4*Fs) );
                end
                
                S(i)=nirs.core.ImageStats;
                S(i).description = data(i).description;
                S(i).demographics   = data(i).demographics;
                S(i).mesh=data(i).mesh.mesh;
                S(i).dfe  = stats.dfe(1);
                
                nVox=size(data(i).mesh.mesh.nodes);
                if(nVox==0);
                    nVox=max(data(i).mesh.link.vertex);
                end
                
                Vall=sparse(nVox,size(data(i).data,2));
                Vall(data(i).mesh.link.vertex,:)=data(i).projectors;
                
                ncond = length(names);
                Vall=repmat(Vall,ncond,1);
               
                nchan=size(Vall,2);
                S(i).beta=vec(stats.beta(1:ncond,:)*Vall');
                covb = sparse(ncond*size(Vall,2),ncond*nchan);
                for id=1:ncond
                    for id2=1:ncond
                        covb((id-1)*nchan+[1:nchan],(id2-1)*nchan+[1:nchan])=sparse(W*diag(squeeze(stats.covb(id,id2,:)))*W');
                    end
                end
                [Uu,Su,Vu]=nirs.math.mysvd(covb);
                S(i).covb_chol = Vall*Uu*sqrt(Su);  % Note- SE = sqrt(sum(G.covb_chol.^2,2))
                
                flds=unique(data(i).mesh.link.type);
                tbl=[];
                for id=1:length(names)
                    for j=1:length(flds)
                        tbl=[tbl; table(arrayfun(@(x){x},[1:nVox]'),...
                            repmat({flds{j}},nVox,1),...
                            repmat({names{id}},nVox,1),...
                            'VariableNames',{'VoxID','type','cond'})];
                    end
                end
                S(i).variables=tbl;
                
                
     
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

