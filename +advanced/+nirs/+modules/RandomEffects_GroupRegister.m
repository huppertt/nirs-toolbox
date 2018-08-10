classdef RandomEffects_GroupRegister < nirs.modules.AbstractGLM
    %% RandomEffects_GroupRegister
    
    % Options:
    %     basis       - a Dictionary object containing temporal bases using stim name as key
    %     verbose     - flag to display progress
    %     trend_func  - a function that takes in a time vector and returns trend regressors
    %     type        - {OLS, NIRS-SPM, or [AR-IRLS]}
    % Example:
    %     j = nirs.modules.GLM();
    %     j.type = 'AR-IRLS';
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
        type;
        options;
        sortfield;
        registration_type;
    end
    methods
        function obj = RandomEffects_GroupRegister( prevJob )
            if nargin > 0, obj.prevJob = prevJob; end
            
            obj.name = 'RandomEffects_GroupRegister';
            obj.basis('default') = nirs.design.basis.Canonical();
            obj.type='AR-IRLS';
            obj.options=[];
            obj.sortfield='subject';
            obj.registration_type='rigid';
        end
        
        function obj = set.type(obj,type)
            validtypes={'OLS','NIRS-SPM','AR-IRLS','MV-GLM','Nonlinear'};
            if(~ismember(type,validtypes))
                disp('type must be one of : ')
                disp(strvcat(validtypes));
                return;
            else
                obj.type=type;
            end
            % use the call functions to evoke any special messages or
            % conditions (e.g. NIRS-SPM checks for code in path);
            switch(obj.type)
                case('OLS')
                    j=nirs.modules.OLS();
                case('AR-IRLS');
                    j=nirs.modules.AR_IRLS();
                case('NIRS-SPM')
                    j=nirs.modules.NIRS_SPM_GLM();
                case('MV-GLM')
                    disp(['Inputs expected to be optical density']);
                case('Nonlinear')
                    j=nirs.modules.nonlin_GLM();
                    obj.basis=j.basis;
                otherwise
                    error('type not recognized');
            end
            obj.options=[];
            pj=obj.prevJob;
            flds=fields(j);
            lst=find(ismember(flds,fields(obj)));
            for i=1:length(lst)
                obj.(flds{lst(i)})=j.(flds{lst(i)});
            end
            lst=find(~ismember(flds,fields(obj)));
            for i=1:length(lst)
                obj.options=setfield(obj.options,flds{lst(i)},j.(flds{lst(i)}));
            end
            obj.prevJob=pj;
        end
        
        function data = runThis( obj, data )
            switch(obj.type)
                case('OLS')
                    j=nirs.modules.OLS();
                case('AR-IRLS');
                    j=nirs.modules.AR_IRLS();
                case('NIRS-SPM')
                    j=nirs.modules.NIRS_SPM_GLM();
                case('MV-GLM')
                    j=nirs.modules.MultiVarGLM();
                case('Nonlinear')
                    j=nirs.modules.nonlin_GLM();
                otherwise
                    error('type not recognized');
            end
            
            j.basis=obj.basis;
            j.verbose=obj.verbose;
            j.trend_func=obj.trend_func;
            j.goforit=obj.goforit;
            
            if(~isempty(obj.options))
                flds=fields(obj.options);
                for i=1:length(flds)
                    j.(flds{i})=obj.options.(flds{i});
                end
            end
            
            S=j.run(data);
            
            j=nirs.modules.SubjLevelStats;
            j.sortfield=obj.sortfield;
            S=j.run(S);
            
            j=nirs.modules.MixedEffects;
            j.formula='beta ~ -1 + cond';
            SOrg=S;
            for iter=1:2
                disp(['Iteration ' num2str(iter) ' of 2']);
                for i=1:length(S)
                    G=j.run(S([1:i-1 i+1:end]));
                    [S(i) tform(i,iter)]=nirs.classify.util.warp_stats(G,S(i),obj.registration_type);
                end
            end
            for iter=1:2
                for i=1:length(S)
                    lst=find(ismember(nirs.createDemographicsTable(data).(obj.sortfield),...
                        nirs.createDemographicsTable(S(i)).(obj.sortfield)));
                    for j=1:length(lst)
                        data(lst(j))=nirs.classify.util.apply_warp(data(i),tform(i,iter));
                    end
                end
            end
            
            
        end
        
    end
    
end

