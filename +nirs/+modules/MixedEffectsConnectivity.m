classdef MixedEffectsConnectivity < nirs.modules.AbstractModule
    %% MixedEffect - Performs group level mixed effects analysis.
    %
    % Options:
    %     formula     - string specifiying regression formula (see Wilkinson notation)
    %     dummyCoding - dummyCoding format for categorical variables (full, reference, effects)
    %     centerVars  - (true or false) flag for whether or not to center numerical variables
    %
    % Example Formula:
    %     % this will calculate the group average for each condition
    %     j = nirs.modules.MixedEffectsConnectivity();
    %     j.formula = 'R ~ -1 + group:cond + (1|subject)';
    %     j.dummyCoding = 'full';
    
    properties
        formula = 'R ~ -1 + cond';
        dummyCoding = 'full';
        centerVars = true;
        robust = false;
        include_diagnostics=false;
    end
    
    methods
        function obj = MixedEffectsConnectivity( prevJob )
            obj.name = 'Mixed Effects Model for Connectivity';
            if nargin > 0
                obj.prevJob = prevJob;
            end
        end
        
        function G = runThis( obj, S )
            
            % demographics info
            demo = nirs.createDemographicsTable( S );
            
            if(isempty(demo))
                for idx=1:length(S)
                    if(iscell(S(idx).demographics))
                        %hyperscanning
                        S(idx).demographics=Dictionary;
                    end
                    
                    S(idx).demographics('fileIdx')=idx;
                end
                demo = nirs.createDemographicsTable( S );
            end
            
            % center numeric variables
            if obj.centerVars
                n = demo.Properties.VariableNames;
                for i = 1:length(n)
                    if all( isnumeric( demo.(n{i}) ) )
                        demo.(n{i}) = demo.(n{i}) - nanmean( demo.(n{i}) );
                    end
                end
            end
            
            % Let's do this per channel for now
            connections=[];
            for i=1:length(S)
                connections=[connections; S(i).probe.connections];
            end
            connections=unique(connections);
            n=height(connections);
            
            if(isa(S(1),'nirs.core.sFCStats'))
                fld='Z';
            else
                fld='F';
            end
            
            D=NaN(length(S),n);
            cond={};
            cnt=1;
            demoall=demo(1,:);
            for i=1:length(S)
                [~,idA]=ismember(S(i).probe.connections,connections);
                tbl=S(i).table;
                for cIdx=1:length(S(i).conditions)
                    lst=find(ismember(tbl.condition,S(i).conditions{cIdx}));
                    D(cnt,idA)=real(tbl.(fld)(lst));
                    cond{cnt,1}=S(i).conditions{cIdx};
                    bind = strfind(cond{cnt,1},'◄');
                    if ~isempty(bind)
                        cond{cnt,1} = cond{cnt,1}(1:bind-2);
                    end
                    demoall(cnt,:)=demo(i,:);
                    cnt=cnt+1;
                end
            end
            % 
            % mask = triu(true(sqrt(n)),1);
            % for i=1:length(S)
            % 
            %     for cIdx=1:length(S(i).conditions)
            %         D(cnt,:)=real(reshape(S(i).(fld)(:,:,cIdx),[],1));
            %         cond{cnt,1}=S(i).conditions{cIdx};
            %         bind = strfind(cond{cnt,1},'◄');
            %         if ~isempty(bind)
            %             cond{cnt,1} = cond{cnt,1}(1:bind-2);
            %         end
            %         demoall(cnt,:)=demo(i,:);
            %         cnt=cnt+1;
            %     end
            % end
            D(D==Inf)=1/eps(1);
            D(D==-Inf)=-1/eps(1);
            
            NoInter=[];
            if(~isempty(strfind(obj.formula,'{')))
                % the formula has variables of no interest
                lstt=sort([strfind(obj.formula,'{') strfind(obj.formula,'}')]);
                cnt=1;
                for ii=1:2:length(lstt)
                    NoInter{cnt}=obj.formula([lstt(ii)+1:lstt(ii+1)-1]);
                    cnt=cnt+1;
                end
                obj.formula=strrep(obj.formula,'{',' ');
                obj.formula=strrep(obj.formula,'}',' ');
            end

            formula=obj.formula;
            formula=['corr ' formula(strfind(formula,'~'):end)];
            nRE=length(strfind(obj.formula,'|'));
            
            corr=rand(size(D(:,1)));
            vars=[table(corr,cond) demoall];
            warning('off','stats:classreg:regr:lmeutils:StandardLinearMixedModel:Message_PerfectFit');
            
            lm = fitlme(vars,formula, 'dummyVarCoding',obj.dummyCoding,...
                'FitMethod', 'ML');
            X = lm.designMatrix('Fixed');
            Z = lm.designMatrix('Random');
            
            % Fit the model
            [Coef,~,StdErr] = nirs.math.fitlme(X,D,Z,obj.robust,false,false);
            
            if(obj.include_diagnostics)
                if(obj.robust)
                    robustflag='on';
                else
                    robustflag='off';
                end
                varnames=lm.CoefficientNames'; %VariableNames(find(lm.VariableInfo.InModel))';
                varnames{end+1}='corr';
                for idx=1:size(D,2)
                    warning('off','stats:statrobustfit:IterationLimit');
                    models{idx} = fitlm(X,D(:,idx),'linear','RobustOpts',robustflag,'Intercept',false,'VarNames',varnames);
                end
                assignin('base','ME_Conn_models',models);
                disp('Diagnostics variable created in workspace named: ME_Conn_models')
            end
             

            % Get results into correct layout
            Coef = Coef(:); %permute(reshape(Coef,[size(Coef,1) sqrt(n) sqrt(n)]),[2 3 1]);
            StdErr = permute(StdErr,[3 1 2]); % permute(reshape(StdErr,[size(StdErr,1) size(StdErr,2) sqrt(n) sqrt(n)]),[3 4 1 2]);
            
            CoefficientNames=lm.CoefficientNames;
            
            %Now sort back out
            G = nirs.core.sFCStats();
            G.type=S(1).type;
            
            % Strip out coefficient name prefixes, e.g. "cond_Puzzle:group_Control" -> "Puzzle:Control"
            PredictorNames=lm.PredictorNames;
            for i=1:length(CoefficientNames)
                CoeffParts = strsplit(CoefficientNames{i},':');
                for j = 1:length(CoeffParts)
                    for k = 1:length(PredictorNames)
                        if ~isempty(strfind( CoeffParts{j} , PredictorNames{k} ))
                            CoeffParts{j} = strrep( CoeffParts{j} , [PredictorNames{k} '_'] , '' );
                        end
                    end
                end
                CoefficientNames{i} = strjoin( CoeffParts , ':' );
            end
            
            %G.conditions = CoefficientNames;
            G.description = 'Group Level Connectivity';
            G.probe=S(1).probe;   

            Cnames=reshape(repmat(CoefficientNames',height(connections),1),[],1);
            connections=repmat(connections,length(CoefficientNames),1);
            connections.type=Cnames;
            G.probe.connections=connections;
            %[n,m]=size(S(1).R);
            G.R=tanh(Coef);
            G.ZstdErr = StdErr;
            G.dfe=lm.DFE;

            

            
            %G.variables.model=models;

             %Remove variables of no interest
            if(~isempty(NoInter))
                
                PredictorNames=lm.PredictorNames;
                tmp=vars;
                for idd=1:length(PredictorNames);
                    if(~isnumeric(tmp.(PredictorNames{idd})))
                        upred=unique(tmp.(PredictorNames{idd}));
                        NoInter=repmat(NoInter,length(upred)*2,1);
                        for ii=1:length(upred)
                            for jj=1:size(NoInter,2)
                                NoInter{ii,jj}=strrep(NoInter{ii,jj},PredictorNames{idd},upred{ii});
                            end
                        end
                        for ii=1:length(upred)
                            for jj=1:size(NoInter,2)
                                NoInter{ii+length(upred),jj}=strrep(NoInter{ii+length(upred),jj},PredictorNames{idd},[PredictorNames{idd} '_' upred{ii}]);
                            end
                        end
                        NoInter=unique(NoInter(:));
                    end
                end
                
                cnames = G.conditions(:);

                cnames=unique(cnames);
                remove={};
                for i=1:length(NoInter); 
                    ss=strsplit(NoInter{i},':'); 
                    for jj=1:length(cnames)
                        flag=true;
                        for ii=1:length(ss)
                            flag=flag & contains(cnames{jj},ss{ii});
                        end
                        if(flag)
                            remove{end+1}=cnames{jj};
                            disp(['Removing condition of no interest: ' remove{end}]);
                        end
                    end
                end;
                if(~isempty(remove))
                    job=nirs.modules.DiscardStims;
                    job.listOfStims=remove;
                    G=job.run(G);

                end
            end
            
        end
    end
    
end