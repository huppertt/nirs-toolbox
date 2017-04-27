classdef GLM < nirs.modules.AbstractGLM
%% GLM- this is a wrapper for the OLS, AR-IRLS, and SPM-NIRS regression programs%
%
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
    end
    methods
        function obj = GLM( prevJob )
            if nargin > 0, obj.prevJob = prevJob; end
            
            obj.name = 'GLM model';
            obj.basis('default') = nirs.design.basis.Canonical();
            obj.type='AR-IRLS';
            obj.options=[];
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
        
        function S = runThis( obj, data )
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
                    
                    
            
        end
        
    end
    
end

