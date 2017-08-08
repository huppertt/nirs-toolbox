classdef (AllowedSubclasses = {?GeneralizedLinearModel, ?LinearModel}) TermsRegression < classreg.regr.ParametricRegression
%TermsRegression Fitted predictive regression model based on terms.
%   TermsRegression is an abstract class representing a fitted regression
%   model for predicting a response as a linear function terms computed
%   from predictor variables. You cannot create instances of this class
%   directly.  You must create a derived class by calling the fit method of
%   a derived class such as LinearModel, GeneralizedLinearModel, or
%   NonLinearModel.
%
%   See also LinearModel, GeneralizedLinearModel, NonLinearModel.

%   Copyright 2011-2014 The MathWorks, Inc.

properties(GetAccess='protected',SetAccess='protected')
        Leverage = [];
        DummyVarCoding = 'reference';
end

properties(GetAccess='public',SetAccess='protected')
%Steps - Structure of stepwise regression results.
%   If the model was fit using stepwise regression, the Steps property
%   is a structure containing information about that fit. If the model was
%   not fit using stepwise regression, this property is empty.
%
%   The Steps structure contains the following fields:
%      Start    Formula representing the initial fit
%      Lower    Formula representing terms that must remain in the fit
%      Upper    Formula representing terms available to the fit
%      Criterion Criterion used for the fit, such as 'SSE'
%      PEnter   P-value below which terms are added to the model
%      PRemove  P-value above which terms are removed from the model
%      History  Table representing the steps taken in the fit
%
%   The History table has one row for each step counting the initial fit,
%   and the following variables (columns):
%      Action   Action taken: Start, Add, or Remove
%      TermName Term moved in that step
%      Terms    Matrix of all terms include in the fit at that step
%      DF       Degrees of freedom for the fit in that step
%      delDF    Change in degrees of freedom at that step
%      Deviance Deviance (residual sum of squares) at that step
%      FStat    F statistic for the term added or removed at that step
%      PValue   P-value for the F statistic
%
%   See also TermsRegression.
        Steps = [];
    end
    properties(GetAccess='protected',SetAccess='protected')
        Design = []; % not reduced with Subset
        CoefTerm = []; % what term is each coefficient (column of Design) part of?
    end
    properties(Dependent,GetAccess='protected',SetAccess='protected')
        % "Working" values - created and saved during fit, but subject to
        % being cleared out and recreated when needed
        design_r = []; % reduced, i.e. Design(Subset,:)
    end
    properties(Constant,Access='protected')
        DummyVarCodings = {'full' 'reference' 'referenceLast' 'effects' 'difference' 'backwardDifference'};
    end
    
    methods(Access='protected')
        function H = get_HatMatrix(model)
            if hasData(model)
                % Let Xw = sqrt(W)*X.
                % yfit = X * inv(X'*W*X) * X'*W*y 
                %      = sqrt(W^-1)* [Xw*inv(Xw'*Xw)*Xw] *sqrt(W)*y
                %      = H*y
                % Compute [...] using a QR decomposition on Xw.
                w_r = get_CombinedWeights_r(model);
                sw = sqrt(w_r);
                X_r = model.design_r;
                Xw_r = bsxfun(@times,X_r,sw);
                [Qw,~,~] = qr(Xw_r,0); % do the same pivoting as in the real fit
                rank = model.NumEstimatedCoefficients;
                Qw = Qw(:,1:rank);
                T = Qw*Qw';
                H = zeros(model.NumObservations);
                subset = model.ObservationInfo.Subset;
                H1 = bsxfun(@times,bsxfun(@times,1./sw,T),sw');
                H1(sw<=0,:) = 0; % prefer 0 to NaN for pt with 0 weight
                H(subset,subset) = H1;
            else
                H = [];
            end
        end
        function d = get_CooksDistance(model)
            if hasData(model)
                w = get_CombinedWeights_r(model,false);
                r = model.Residuals.Raw;
                h = model.Leverage;
                d = w .* abs(r).^2 .* (h./(1-h).^2)./(model.NumEstimatedCoefficients*varianceParam(model));
            else
                d = [];
            end
        end
        function checkDesignRank(model)
            if model.NumEstimatedCoefficients<model.NumCoefficients
                warning(message('stats:LinearModel:RankDefDesignMat'));
            end
        end
        function w = get_CombinedWeights_r(model,reduce)
            % classes with no special (robust, glm) weights can use this
            w = model.ObservationInfo.Weights;
            if nargin<2 || reduce
                subset = model.ObservationInfo.Subset;
                w = w(subset);
            end
        end
    end
        
    methods % get/set methods
        function design_r = get.design_r(model)
            if isempty(model.WorkingValues)
                design_r = create_design_r(model);
            else
                design_r = model.WorkingValues.design_r;
            end
        end
    end % get/set methods
    
    methods(Access='public')
        function model = addTerms(model,terms)
%addTerms Add terms to a regression model.
%   M2 = addTerms(M1,TERMS) adds the terms specified by TERMS to the
%   regression model M1 and returns the new fitted model as M2. TERMS can
%   be a text string representing one or more terms (as used on the
%   right-hand-side of a formula supplied to the FIT method, or a terms
%   matrix.
%
%   Example:
%       % Fit model to car data; add interaction term
%       load carsmall
%       d = dataset(MPG,Weight);
%       d.Year = ordinal(Model_Year);
%       lm1 = fitlm(d,'MPG ~ Year + Weight + Weight^2')
%       lm2 = addTerms(lm1,'Year*Weight')
%
%   See also LinearModel, GeneralizedLinearModel, removeTerms.

            model.Formula = addTerms(model.Formula,terms);
            model = removeCategoricalPowers(model);
            model = doFit(model);
            checkDesignRank(model)
        end
        function model = removeTerms(model,terms)
%removeTerms Remove terms from a regression model.
%   M2 = removeTerms(M1,TERMS) removes the terms specified by TERMS from
%   the regression model M1 and returns the new fitted model as M2. TERMS
%   can be a text string representing one or more terms (as used on the
%   right-hand-side of a formula supplied to the FIT method, or a terms
%   matrix.
%
%   Example:
%       % Fit model to car data; remove squared term
%       load carsmall
%       d = dataset(MPG,Weight);
%       d.Year = ordinal(Model_Year);
%       lm1 = fitlm(d,'MPG ~ Year + Weight + Weight^2')
%       lm2 = removeTerms(lm1,'Weight:Weight')
%
%   See also LinearModel, GeneralizedLinearModel, addTerms.

            model.Formula = removeTerms(model.Formula,terms);
            model = doFit(model);
            checkDesignRank(model)
        end
        
        function model = step(model,varargin)
            paramNames = {   'Lower'        'Upper'  'Criterion' 'PEnter' 'PRemove' 'NSteps' 'Verbose'};
            paramDflts = {'constant' 'interactions'        'SSE'       []        []        1         1};
            
            wasempty = isempty(model.Steps);
            if ~wasempty
                paramDflts{1} = model.Steps.Lower.Terms;
                paramDflts{2} = model.Steps.Upper.Terms;
                paramDflts{3} = model.Steps.Criterion;
                paramDflts{4} = model.Steps.PEnter;
                paramDflts{5} = model.Steps.PRemove;
            end
            
            [lower,upper,crit,penter,premove,nsteps,verbose] = ...
                internal.stats.parseArgs(paramNames, paramDflts, varargin{:});

            if ~isscalar(verbose) || ~ismember(verbose,0:2)
                error(message('stats:LinearModel:BadVerbose'));
            end
            
            start = model.Formula;
            
            if ~isa(lower,'LinearFormula')
                lower = classreg.regr.LinearFormula(lower,start.VariableNames,start.ResponseName,start.HasIntercept,start.Link);
            end
            if ~isa(upper,'LinearFormula')
%                 if classreg.regr.LinearFormula.isModelAlias(upper)
%                     upper = {upper,start.PredictorNames};
%                 end
                upper = classreg.regr.LinearFormula(upper,start.VariableNames,start.ResponseName,start.HasIntercept,start.Link);
            end
            
%             [penter,premove] = classreg.regr.TermsRegression.getDefaultThresholds(crit,penter,premove);
%
            model.Steps.Start = start;
            model.Steps.Lower = lower;
            model.Steps.Upper = upper;
            model.Steps.Criterion = crit;
            model.Steps.PEnter = penter;
            model.Steps.PRemove = premove;
            if wasempty
                model.Steps.History = [];
            end

            model = stepwiseFitter(model,nsteps,verbose);
        end
        
        function [p,t,r] = coefTest(model,H,c)
%coefTest Linear hypothesis test on coefficients.
%   P = coefTest(M) computes the p-value for an F test that all
%   coefficient estimates in the regression model M except the intercept
%   are zero.
%
%   P = coefTest(M,H), with H a numeric matrix having one column for each
%   coefficient, performs an F test that H*B=0, where B represents the
%   coefficient vector.
%
%   P = coefTest(M,H,C) accepts a vector C having one value for each row
%   of H, and it performs an F test that H*B=C.
%
%   [P,F,R] = coefTest(...) also returns the F-statistic F and the rank R
%   of the matrix H. The F statistic has R degrees of freedom in the
%   numerator and M.DFE degrees of freedom in the denominator.
%
%   Example:
%       % Test the significance of the Weight^2 coefficient, and note that
%       % the p-value is the same as in the coefficients display.
%       load carsmall
%       d = dataset(MPG,Weight);
%       d.Year = ordinal(Model_Year);
%       lm = fitlm(d,'MPG ~ Year + Weight + Weight^2')
%       p = coefTest(lm,[0 0 0 0 1])
%
%       % Test the significance of both coefficients for the categorical
%       % predictor Year, and note that the p-value is the same as in the
%       % anova table display.
%       anova(lm)
%       p = coefTest(lm,[0 0 1 0 0;0 0 0 1 0])
%
%   See also LinearModel, GeneralizedLinearModel, NonLinearModel, linhyptest.

            nc = model.NumCoefficients;
            if nargin<2
                % Default is to test all terms except the constant
                H = eye(nc);
                constTerm = find(all(model.Formula.Terms==0,2));
                if isscalar(constTerm)
                    coefRow = model.CoefTerm==constTerm;
                    H(coefRow,:) = [];
                end
            end
            if nargin < 3
                c = zeros(size(H,1),1);
            end
            [p,t,r] = coefTest@classreg.regr.ParametricRegression(model,H,c);
        end
    end % public
    
    methods(Access='protected')
        function model = TermsRegression()
            model.PredictorTypes = 'mixed';
            if nargin == 0 % special case
                model.Formula = classreg.regr.LinearFormula;
                return
            end
        end
        
        function D = get_diagnostics(model,type)
            if nargin<2 % return all diagnostics in a table
                CooksDistance = get_diagnostics(model,'cooksdistance');
                HatMatrix = get_diagnostics(model,'hatmatrix');
                Leverage = model.Leverage;
                D = table(Leverage,CooksDistance,HatMatrix,...
                            'RowNames',model.ObservationNames);
            else        % return a single diagnostic
                subset = model.ObservationInfo.Subset;
                switch(lower(type))
                    case 'leverage'
                        D = model.Leverage;
                        D(~subset,:) = 0;                  
                    case 'hatmatrix'
                        try
                            D = get_HatMatrix(model);
                        catch ME
                            warning(message('stats:LinearModel:HatMatrixError', ...
                                ME.message));
                            D = zeros(length(subset),0);
                        end
                        D(~subset,:) = 0;
                    case 'cooksdistance'
                        D = get_CooksDistance(model);
                        D(~subset,:) = NaN;
                    otherwise
                        error(message('stats:LinearModel:UnrecognizedDiagnostic', type));
                end
            end
        end
         
        % --------------------------------------------------------------------
        function [design,coefTerm,coefNames] = designMatrix(model,X,isPredictorsOnly,terms)
            % Translate a data matrix into a design matrix.
            if nargin<4
                terms = model.Formula.Terms;
            end
            if isa(X,'dataset')
                X = dataset2table(X);
            end
                
            if isa(X,'table')
                % Given a dataset/table, make sure it has all of the predictor
                % variables.  It can also have any unused variables from
                % the model fit.
                tf = ismember(model.PredictorNames,X.Properties.VariableNames);
                if ~all(tf)
                    error(message('stats:classreg:regr:TermsRegression:MissingVariable'));
                end
                [tf,varLocs] = ismember(X.Properties.VariableNames,model.VariableNames);
                if ~all(tf)
                    % Remove extra variables from X so that we don't have to
                    % construct Terms, IsCategorical, Range, and
                    % DummyVarCoding values.
                    X = X(:,tf);
                end
                varLocs = varLocs(tf);
                [design,~,~,coefTerm,coefNames] ...
                    = classreg.regr.modelutils.designmatrix(X,'Model',terms(:,varLocs), ...
                    'DummyVarCoding',model.DummyVarCoding, ...
                    'CategoricalVars',model.VariableInfo.IsCategorical(varLocs), ...
                    'CategoricalLevels',model.VariableInfo.Range(varLocs));
            else
                nvars = model.NumVariables;
                if nargin > 2 && ~isempty(isPredictorsOnly) && isPredictorsOnly
                    % Given a matrix, assume the columns correspond in order
                    % to the predictor variables
                    npreds = model.NumPredictors;
                    varLocs = model.PredLocs;
                else
                    % Given a matrix, assume the columns correspond in order
                    % to the non-response variables
                    npreds = nvars-1;
                    varLocs = true(1,nvars);
                    varLocs(model.RespLoc) = false;
                end
                if size(X,2) ~= npreds
                    error(message('stats:classreg:regr:TermsRegression:WrongXColumns', nvars - 1));
                end
                [design,~,~,coefTerm,coefNames] ...
                    = classreg.regr.modelutils.designmatrix(X,'Model',terms(:,varLocs), ...
                    'CategoricalVars',model.VariableInfo.IsCategorical(varLocs), ...
                    'CategoricalLevels',model.VariableInfo.Range(varLocs), ...
                    'DummyVarCoding',model.DummyVarCoding, ...
                    'VarNames',model.Formula.VariableNames(varLocs));
            end
        end
        
        % --------------------------------------------------------------------
        function model = assignData(model,X,y,w,asCat,dummyCoding,varNames,excl)
            model = assignData@classreg.regr.ParametricRegression(model,X,y,w,asCat,varNames,excl);
            
            tf = internal.stats.isString(dummyCoding);
            if ~tf || ~ismember(dummyCoding,classreg.regr.TermsRegression.DummyVarCodings)
                error(message('stats:classreg:regr:TermsRegression:BadCodingValue', internal.stats.listStrings( classreg.regr.TermsRegression.DummyVarCodings )));
            end
            model.DummyVarCoding = dummyCoding;
        end
        
        % --------------------------------------------------------------------
        function model = selectVariables(model)
            f = model.Formula;
            [~,model.PredLocs] = ismember(f.PredictorNames,f.VariableNames);
            [~,model.RespLoc] = ismember(f.ResponseName,f.VariableNames);
            model = selectVariables@classreg.regr.ParametricRegression(model);
        end
        
        % --------------------------------------------------------------------
        function model = postFit(model)
            model = postFit@classreg.regr.ParametricRegression(model);
            wts = get_CombinedWeights_r(model);
            Xw_r = bsxfun(@times,model.design_r,sqrt(wts));
            [Qw_r,~,~] = qr(Xw_r,0); % do the same pivoting as in the real fit
            h = zeros(size(model.ObservationInfo,1),1);
            h(model.ObservationInfo.Subset) = sum(abs(Qw_r).^2,2);
            model.Leverage = h;
        end
        
        % --------------------------------------------------------------------
        function design_r = create_design_r(model)
            design_r = model.Design(model.ObservationInfo.Subset,:);
        end
        
        % --------------------------------------------------------------------
        function model = removeCategoricalPowers(model,silent)
            if nargin<2
                silent = false;
            end
            f = model.Formula;
            isCat = model.VariableInfo.IsCategorical;
            model.Formula = removeCategoricalPowers(f,isCat,silent);
        end
        
        % --------------------------------------------------------------------
        function tf = hasConstantModelNested(model)
            % a model with an explicit intercept, or with a categorical
            % variable with full dummy variable encoding, contains the
            % constant model
            tf = model.Formula.HasIntercept;
            if ~tf && hasData(model)
                catVarsInModel = model.VariableInfo.InModel & model.VariableInfo.IsCategorical;
                tf = strcmpi('full',model.DummyVarCoding) && any(catVarsInModel);
            end
        end
        
        function fit = stepwiseFitter(start,nsteps,verbose)
            lowerBound = start.Steps.Lower.Terms;
            upperBound = start.Steps.Upper.Terms;
            crit = start.Steps.Criterion;
            [addTest,addThreshold,removeTest,removeThreshold,reportedNames,testName] = ...
                start.getStepwiseTests(crit,start.Steps);
            
            terms = start.Formula.Terms;
            varNames = start.Formula.VariableNames;
            [~,~,startReportedVals] = addTest(start,[]);
            if isempty(start.Steps.History)
                history = table(nominal('Start',[],{'Start' 'Add' 'Remove'}), ...
                    {start.Formula.LinearPredictor}, ...
                    {terms}, start.NumEstimatedCoefficients,NaN,startReportedVals{:}, ...
                    'VariableNames',[{'Action', 'TermName' 'Terms' 'DF' 'delDF'} reportedNames]);
            else
                history = start.Steps.History;
            end
            fit = start;
            justAdded = [];
            justRemoved = [];
            
            while nsteps > 0 % nsteps might be inf
                nsteps = nsteps - 1;
                changed = false;
                
                % Try to add a term to the model
                candidates = find(candidatesToAdd(terms,upperBound,justRemoved));
                bestTestVal = Inf;
                X = getData(fit);
                [Qdesign,~] = qr(fit.design_r,0);
                inclRows = fit.ObservationInfo.Subset;
                for j = 1:length(candidates)
                    newtermj = upperBound(candidates(j),:);
                    [termsj,ord] = sortTerms([terms; newtermj]);
                    [~,locj] = max(ord); % where did upperBound(candidates(j),:) end up?
                    newxterms = designMatrix(fit,X,[],newtermj);

                    if redundantTerm(Qdesign,newxterms,inclRows)
                        testValj = Inf;
                    else
                        fitj = reFit(fit,termsj);
                        % Remember the best model in this go-round
                        [testValj,testValjReported,reportedValsj] = addTest(fitj,fit);
                        if verbose > 1
                            tn = classreg.regr.modelutils.terms2names(termsj(locj,:),varNames); tn = tn{:};
                            fprintf('   %s',getString(message('stats:classreg:regr:TermsRegression:display_ForAddingIs',testName,tn,num2str(testValjReported))));
                        end
                    end
                    if testValj < bestTestVal
                        bestj = j;
                        bestTerms = termsj;
                        bestLoc = locj; % location in bestTerms of best term so far
                        bestFit = fitj;
                        bestTestVal = testValj;
                        bestReportedVals = reportedValsj;
                    end
                end
                if verbose > 1 && isempty(candidates)
                    fprintf('   %s',getString(message('stats:classreg:regr:TermsRegression:display_NoCandidateTermsToAdd')));
                end
                % Update the model if adding the term improves it sufficiently
                if bestTestVal < addThreshold
                    addedTermName = classreg.regr.modelutils.terms2names(bestTerms(bestLoc,:),varNames);
                    delDF = fit.DFE - bestFit.DFE;
                    history(end+1,:) = table(nominal('Add'),addedTermName,{bestTerms},bestFit.NumEstimatedCoefficients,delDF,bestReportedVals{:});
                    if verbose
                        allVals = strcat(reportedNames(:),{' = '},strjust(num2str(cat(1,bestReportedVals{:})),'left'));
                        displayString = sprintf(', %s',allVals{:});
                        fprintf('%d. Adding %s%s\n',size(history,1)-1,history.TermName{end},displayString);
                    end
                    terms = bestTerms;
                    fit = bestFit;
                    changed = true;
                    justAdded = bestLoc;
                    justRemoved = [];
                    
                else
                    if verbose > 2
                        fprintf('   %s',getString(message('stats:classreg:regr:TermsRegression:display_NoTermsToAddSmallestPValue',num2str(bestTestVal))));
                    end
                    % Try to remove a term from the model
                    candidates = find(candidatesToRemove(terms,lowerBound,justAdded));
                    bestTestVal = -Inf;
                    terminfo = getTermInfo(fit);
                    designterms = terminfo.designTerms;
                    for j = 1:length(candidates)
                        termsj = terms;
                        termsj(candidates(j),:) = [];
                        dtj = (designterms == candidates(j));
                        newxterms = fit.design_r(:,dtj);
                        [qrd,~] = qr(fit.design_r(:,~dtj),0);
                        fitj = reFit(fit,termsj);
                        % Remember the best model in this go-round
                        [testValj,testValjReported,reportedValsj] = removeTest(fitj,fit);
                        if verbose > 1
                            tn = classreg.regr.modelutils.terms2names(terms(candidates(j),:),varNames); tn = tn{:};
                            fprintf('   %s',getString(message('stats:classreg:regr:TermsRegression:display_ForRemovingIs',testName,tn,num2str(testValjReported))));
                        end

                        % Sometimes the following is helpful to get rid of
                        % terms that cause a singular model
%                         if redundantTerm(qrd,newxterms)
%                             testValj = NaN;
%                         end
                        if isnan(testValj) || testValj>bestTestVal
                            % We're taking NaN as a sign of a term that is
                            % somehow bad and ought to be removed
                            bestj = j;
                            bestTerms = termsj;
                            bestFit = fitj;
                            bestTestVal = testValj;
                            bestReportedVals = reportedValsj;
                            if isnan(testValj)  % ready to remove this now
                                break
                            end
                        end
                    end
                    if verbose > 1 && isempty(candidates)
                        fprintf('   %s',getString(message('stats:classreg:regr:TermsRegression:display_NoCandidateTermsToRemove')));
                    end
                    
                    % Update the model if removing the term improves it sufficiently
                    if isnan(bestTestVal) || bestTestVal>removeThreshold
                        removedTerm = terms(candidates(bestj),:);
                        removedTermName = classreg.regr.modelutils.terms2names(removedTerm,varNames);
                        delDF = fit.DFE - bestFit.DFE;
                        history(end+1,:) = table(nominal('Remove'),removedTermName,{bestTerms},bestFit.NumEstimatedCoefficients,delDF,bestReportedVals{:});
                        if verbose
                            allVals = strcat(reportedNames(:),{' = '},strjust(num2str(cat(1,bestReportedVals{:}),'%.5g'),'left'));
                            displayString = sprintf(', %s',allVals{:});
                            fprintf('%s',getString(message('stats:classreg:regr:TermsRegression:display_Removing',size(history,1)-1,history.TermName{end},displayString)));
                        end
                        terms = bestTerms;
                        fit = bestFit;
                        changed = true;
                        justAdded = [];
                        [~,justRemoved] = ismember(removedTerm,upperBound,'rows');
                        if justRemoved==0
                            justRemoved = [];
                        end
                    else
                        if verbose > 2
                            fprintf('%s',getString(message('stats:classreg:regr:TermsRegression:display_DidntRemoveAnything',num2str(bestTestVal))));
                        end
                    end
                end
                if ~changed, break, end
            end
            if verbose && size(history,1)==1
                fprintf(getString(message('stats:classreg:regr:TermsRegression:display_NoTermsToAddToOrRemoveFromInitialModel')));
            end
            fit.Steps = start.Steps;
            % vnames = history.Properties.VarNames(4:end);
            % for i = 1:length(vnames)
            %     vname = vnames{i};
            %     history.(vname) = internal.stats.DoubleTableColumn(history.(vname));
            % end
            fit.Steps.History = history;
        end
        
        function terminfo = getTermInfo(model)
            
            terminfo.terms = model.Formula.Terms;
            terminfo.designTerms = model.CoefTerm;
            terminfo.designMeans = mean(model.Design(model.ObservationInfo.Subset,:),1);
            terminfo.isCatVar = model.VariableInfo.IsCategorical';
            terminfo.numCatLevels = NaN(1,model.NumVariables,1);
            terminfo.numCatLevels(terminfo.isCatVar) = ...
                cellfun(@countlevels,model.VariableInfo.Range(terminfo.isCatVar));
        end
        
    end % protected
    
    methods(Static, Access='protected')
        
        function [penter,premove] = getDefaultThresholds(crit,penter,premove)
            smaller = true; % most criteria are better if smaller

            critString = internal.stats.isString(crit);
            if ~critString && ~isa(crit,'function_handle')
                error(message('stats:classreg:regr:TermsRegression:BadStepwiseCriterion'));
            end
            if (isempty(penter) || isempty(premove)) && ~critString
                error(message('stats:classreg:regr:TermsRegression:MissingThreshold'));
            end

            if critString
                allcrit = {'AIC' 'BIC' 'Rsquared' 'AdjRsquared' 'SSE'};
                crit = internal.stats.getParamVal(crit,allcrit,'''Criterion''');
                
                switch lower(crit)
                    case {'aic' 'bic'}
                        if isempty(penter), penter = 0; end
                        if isempty(premove), premove = 0.01; end
                    case 'rsquared'
                        smaller = false; % better if bigger
                        if isempty(penter), penter = 0.1; end
                        if isempty(premove), premove = 0.05; end
                    case 'adjrsquared'
                        smaller = false; % better if bigger
                        if isempty(penter), penter = 0; end
                        if isempty(premove), premove = -0.05; end
                    case 'sse'
                        if isempty(penter), penter = 0.05; end
                        if isempty(premove), premove = 0.10; end
                end
            end
            
            if smaller && penter>=premove
                error(message('stats:LinearModel:BadSmallerThreshold', sprintf( '%g', penter ), sprintf( '%g', premove ), crit));
            elseif ~smaller && penter<=premove
                error(message('stats:LinearModel:BadLargerThreshold', sprintf( '%g', penter ), sprintf( '%g', premove ), crit));
            end
        end
        
        function [X,y,haveDataset,otherArgs] = handleDataArgs(X,y,varargin) % [X, y | DS], ...
            if isa(X,'dataset');
                X = dataset2table(X);
            end
            
            haveDataset = isa(X,'table');
            if haveDataset
                % Have a dataset/table array, so no y.
                if nargin > 1
                    % Put the second arg in with the rest
                    otherArgs = [{y} varargin];
                else
                    otherArgs = {};
                end
                y = []; % just put anything in y
            elseif nargin < 2
                error(message('stats:classreg:regr:TermsRegression:MissingY'))
            else
                if isrow(X)
                    nx = length(X);
                    if (isvector(y) && numel(y)==nx) || (size(y,1)==nx)
                        X = X';
                    end
                end
                isNumVarX = isnumeric(X) || islogical(X);
                isCatVecX = isa(X,'categorical') && isvector(X);
                if ~(isNumVarX || isCatVecX)
                    error(message('stats:classreg:regr:FitObject:PredictorMatricesRequired'));
                end
                otherArgs = varargin;
            end
        end
        
        function [addTest,addThreshold,removeTest,removeThreshold,reportedNames,testName] = getStepwiseTests(crit,Steps)
            
            addThreshold = Steps.PEnter;
            removeThreshold = Steps.PRemove;
            [addThreshold,removeThreshold] = classreg.regr.TermsRegression.getDefaultThresholds(crit,addThreshold,removeThreshold);
            if internal.stats.isString(crit)
                allcrit = {'AIC' 'BIC' 'Rsquared' 'AdjRsquared' 'SSE'};
                crit = internal.stats.getParamVal(crit,allcrit,'''Criterion''');
                switch lower(crit)
                    case {'aic' 'bic'}
                        addTest = @(proposed,current) ...
                            generic_test(proposed, current, @(fit) get_modelcriterion(fit,crit), 'decreasing');
                        removeTest = @(proposed,current) ...
                            generic_test(proposed, current, @(fit) get_modelcriterion(fit,crit), 'increasing');
                        reportedNames = {upper(crit)};
                        testName = sprintf('%s',getString(message('stats:classreg:regr:TermsRegression:display_ChangeIn',upper(crit))));
                    case {'rsquared' 'adjrsquared'}
                        if strcmpi(crit,'rsquared')
                            rsqType = 'Ordinary';
                        else
                            rsqType = 'Adjusted';
                        end
                        addThreshold = -addThreshold;
                        removeThreshold = -removeThreshold;
                        addTest = @(proposed,current) ...
                            generic_test(proposed, current, @(fit) get_rsquared(fit,rsqType), 'increasing');
                        removeTest = @(proposed,current) ...
                            generic_test(proposed, current, @(fit) get_rsquared(fit,rsqType), 'decreasing');
                        reportedNames = {crit};
                        testName = sprintf('%s',getString(message('stats:classreg:regr:TermsRegression:display_ChangeIn',crit)));
                    case 'sse'
                        addTest = @(proposed,current) f_test(proposed,current,'up');
                        removeTest = @(proposed,current) f_test(current,proposed,'down');
                        reportedNames = {'FStat' 'pValue'};
                        testName = 'pValue';
                end
            elseif isa(crit,'function_handle')
                fun = crit;
                addTest = @(proposed,current) generic_test(proposed, current, fun, 'decreasing');
                removeTest = @(proposed,current) generic_test(proposed, current, fun, 'increasing');
                reportedNames = {'ModelCriterion'};
                testName = getString(message('stats:classreg:regr:TermsRegression:assignment_ChangeInModelCriterion'));
            elseif iscell(crit) && (numel(crit) == 3) && all(cellfun(@(c)isa(c,'function_handle'),crit))
                addTest = crit{1};
                removeTest = crit{2};
                reportedNames = {'ModelCriterion'};
                testName = getString(message('stats:classreg:regr:TermsRegression:assignment_ChangeInModelCriterion'));
            else
                error(message('stats:classreg:regr:TermsRegression:BadStepwiseCriterion'));
            end
        end
        
        function formula = createFormula(supplied,modelDef,X, ...
                predictorVars,responseVar,intercept,link,varNames,haveDataset,clink)
            if nargin<10
                clink = 'identity';
            end
            
            givenTerms = classreg.regr.LinearFormula.isTermsMatrix(modelDef);
            givenAlias = classreg.regr.LinearFormula.isModelAlias(modelDef);
            givenString = ~givenAlias && internal.stats.isString(modelDef);
            
            if ~haveDataset
                if ischar(X)
                    nx = 1;
                else
                    nx = size(X,2);
                end
                [varNames,predictorVars,responseVar] = ...
                        classreg.regr.FitObject.getVarNames(varNames,predictorVars,responseVar,nx);
            end
            
            if isa(modelDef,'classreg.regr.LinearFormula')
                % *** need to compare formula's varspace with X or DS
                formula = modelDef;
                
            elseif givenString || givenTerms || givenAlias
                if haveDataset
                    % VarNames is intended to name the columns of separate X and y.  It's
                    % not needed or accepted for a dataset/table array, because that has names.
                    if supplied.VarNames
                        error(message('stats:classreg:regr:TermsRegression:NoVarNames'));
                    end
                    nvars = size(X,2);
                    varNames = X.Properties.VariableNames;
                else % data given as X,y
                    % ResponseVar is intended to say which variable in a dataset/table array is
                    % the response.  It's not needed or accepted for separate X and y.
                    if supplied.ResponseVar && ~internal.stats.isString(responseVar)
                        error(message('stats:classreg:regr:TermsRegression:BadResponseVar'));
                    end
                    if ischar(X)
                        nvars = 2;  % response plus single predictor
                    else
                        nvars = size(X,2) + 1;
                    end
                    if isempty(varNames)
                        if givenString
                            % If given a formula string and no names, LinearFormula will
                            % assume vars are in alphabetical order, but we have to force the
                            % response to be last.  Since we don't know the response name, we
                            % have to ask LinearFormula first.
                            formula = classreg.regr.LinearFormula(modelDef,[],[],[],link);
                            varNames = formula.VariableNames;
                            varnum = 0; %length(varNames)-1;  % number of x variables
                            while(length(varNames)<nvars)
                                varnum = varnum +1;
                                newname = sprintf('x%d',varnum);
                                if ~ismember(newname,varNames)
                                    varNames{end+1} = newname;
                                end
                            end
                            responseName = formula.ResponseName;
                            respLoc = find(strcmp(responseName,varNames));
                            varNames = varNames([1:(respLoc-1) (respLoc+1):end respLoc]);
                        elseif givenTerms
                            varNames = [];
                        else % givenAlias
                            % LinearFormula will create default names but needs to know
                            % how many variables if only given an alias.
                            varNames = nvars;
                        end
                    else
                        if length(varNames) ~= nvars
                            error(message('stats:classreg:regr:TermsRegression:BadVarNames'));
                        end
                    end
                end
                
                if givenString
                    if ~supplied.Intercept
                        intercept = []; % take the value implied by modelDef
                    end
                    formula = classreg.regr.LinearFormula(modelDef,varNames,'',intercept,link);
                    
                    % If predictor/response vars given, make sure they don't conflict
                    if (supplied.PredictorVars && ~all(ismember(formula.PredictorNames,predictorVars))) || ...
                       (supplied.ResponseVar && ~strcmp(formula.ResponseName,responseVar))
                        error(message('stats:classreg:regr:TermsRegression:NoFormulaVars'));
                    end
                    
                elseif givenTerms
                    if supplied.PredictorVars
                        error(message('stats:classreg:regr:TermsRegression:NoTermsVars'));
                    end
                    if haveDataset
                        if size(modelDef,2) ~= nvars
                            if size(modelDef,2)==1 && size(modelDef,1)==size(X,1)
                                error(message('stats:classreg:regr:TermsRegression:SeparateResponse'));
                            else
                                error(message('stats:classreg:regr:TermsRegression:BadTermsDataset'));
                            end
                        end
                        if ~supplied.ResponseVar
                            responseVar = [];
                        end
                    else % data given as X,y
                        % Already errored if supplied.ResponseVar
                        ncols = size(modelDef,2);
                        if ncols ~= nvars
                            if ncols == nvars-1
                                modelDef(:,nvars) = 0;
                            else
                                error(message('stats:classreg:regr:TermsRegression:BadTermsMatrix'));
                            end
                        end
                        responseVar = nvars;
                    end
                    formula = classreg.regr.LinearFormula(modelDef,varNames,responseVar,intercept,link);
                    
                else % givenAlias
                    if haveDataset
                        if ~supplied.ResponseVar
                            responseVar = [];
                        end
                    elseif ~internal.stats.isString(responseVar) % data given as X,y
                        % Already errored if supplied.ResponseVar
                        responseVar = nvars;
                    end
                    if supplied.PredictorVars
                        modelDef = {modelDef,predictorVars};
                    elseif haveDataset && isempty(predictorVars)
                        predictorVars = getValidPredictors(X,varNames,responseVar);
                        if ~isempty(predictorVars)
                            modelDef = {modelDef,predictorVars};
                        end
                    end
                    
                    formula = classreg.regr.LinearFormula(modelDef,varNames,responseVar,intercept,link);
                end
                
            else
                error(message('stats:classreg:regr:TermsRegression:BadModelDef'));
            end
            
            if ~haveDataset
                if supplied.VarNames && length(formula.VariableNames) ~= nvars
                    error(message('stats:classreg:regr:TermsRegression:BadVarNamesXY'));
                end
            end
            
            % If there is no link in the formula (assumed identity) and nothing
            % in the link argument, use the canonical link
            if isempty(link) && ~isequal(clink,'identity') && isequal(formula.Link,'identity')
                formula.Link = clink;
            end
        end
    end % static protected
    
    methods(Hidden, Static, Access='public')
        function model = loadobj(obj)
            obj = loadobj@classreg.regr.ParametricRegression(obj);
            if isfield(obj.Steps,'History') && isa(obj.Steps.History,'dataset')
                obj.Steps.History = dataset2table(obj.Steps.History);
            end
            model = obj;
        end
    end
end

% ----------------------------------------------------------------------
function [terms,iterms] = sortTerms(terms)
%SORTTERMS Sort linear regression terms
%    SORTEDTERMS = SORTERMS(TERMS) sorts the rows of the linear regression
%    temrs matrix.  Rows of SORTEDTERMS are in nondecreasing term order, that
%    is, SUM(SORTEDTERMS,2) is nondecreasing.  Within each group having the
%    same order, rows of SORTEDTERMS are in nondecreasing order of maximum
%    power, that is, MAX(SORTEDTERMS,2) is nondecreasing within each group.
%    Finally, within each group having the same term order and maximum power,
%    rows of SORTEDTERMS are sorted by the variables that occur in each term.
[nterms,nvars] = size(terms);
iterms = (1:nterms)';
[terms,ord] = unique(terms,'rows'); % remove duplicate terms
iterms = iterms(ord);
[terms,ord] = sortrows(terms,-1:-1:-nvars); % sort by var
iterms = iterms(ord);
[~,ord] = sortrows([sum(terms,2) max(terms,[],2)]);
terms = terms(ord,:); % then sort by term order
iterms = iterms(ord);
end

% ----------------------------------------------------------------------
function candidates = candidatesToAdd(current,upperBound,justRemoved)
% Assumes upperBound is a superset of current
[potentialCandidates,icandidates] = setdiff(upperBound,current,'rows');
candidates = false(size(upperBound,1),1);
candidates(icandidates) = true;

% Don't try to add a term that has just been removed
if ~isempty(justRemoved)
    candidates(justRemoved) = false;
end

% Intercept and linear terms can always be added, but don't add any higher
% order terms for which a lower order subterm has yet to be added to the
% current model.  Note that if a lower order subterm of a potential candidate
% is not present in upperBound (and thus is not in current), that will not
% prevent the candidate from being considered, i.e., there is no assumption
% that current and upperBound are "coherent".
order = sum(potentialCandidates,2);
if max(order) > 1
    [~,ord] = sort(order,1,'ascend');
    not = false(size(potentialCandidates,1),1);
    for i = ord(:)' % check in increasing term order
        subterm = potentialCandidates(i,:);
        termDiffs = bsxfun(@minus,potentialCandidates,subterm);
        j = all(bsxfun(@times,termDiffs,subterm>0)>=0,2);
        j(i) = false;
        not = not | j;
        % Stop checking if no terms are candidates for addition
        if all(not), break; end
    end
    candidates(icandidates(not)) = false;
end
end

% ----------------------------------------------------------------------
function candidates = candidatesToRemove(current,lowerBound,justAdded)
% Assumes lowerBound is a subset of current
[potentialCandidates,icandidates] = setdiff(current,lowerBound,'rows');
candidates = false(size(current,1),1);
candidates(icandidates) = true;

% Don't try to remove a term that has just been added
if ~isempty(justAdded)
    candidates(justAdded) = false;
end

% Don't remove any lower order terms for which a higher order superterm is
% present in the current model.
order = sum(potentialCandidates,2);
[~,ord] = sort(order,1,'descend');
not = false(size(potentialCandidates,1),1);
for i = ord(:)' % check in decreasing term order
    superterm = potentialCandidates(i,:);
    termDiffs = bsxfun(@minus,superterm,potentialCandidates);
    j = all(termDiffs.*(potentialCandidates>0)>=0,2);
    j(i) = false;
    not = not | j;
    % Stop checking if no terms are candidates for removal
    if all(not), break; end
end
candidates(icandidates(not)) = false;
end

% ----------------------------------------------------------------------
function model = reFit(model,terms)
model.Formula.Terms = terms;
model = doFit(model);
end

function [delCrit,critReported,reportedVals] = generic_test(proposed,current,critfun,direction)
% Test function for things like AIC or r-squared where it's a simple
% comparison of the raw model criterion values
newCrit = critfun(proposed);
if isempty(current)
    oldCrit = newCrit;
else
    oldCrit = critfun(current);
end

% Adjust criterion for a "less than is better" test
critReported = newCrit - oldCrit;
if strcmp(direction,'decreasing') % lower critFun value is better
    delCrit = critReported;
else                              % higher critFun value is better
    delCrit = -critReported; 
end
reportedVals = {newCrit};
end

function [p,pReported,reportedVals] = f_test(fit1,fit0,~)
sse1 = fit1.SSE; % the larger model
dfDenom = fit1.DFE;
if isempty(fit0)
    sse0 = sse1;
    dfNumer = 0;
else
    sse0 = fit0.SSE; % the smaller model
    dfNumer = fit0.DFE - fit1.DFE;
end
F = ((sse0-sse1)/dfNumer) / (fit1.SSE/dfDenom);
p = fcdf(1./F,dfDenom,dfNumer); % fpval(F,dfNumer,dfDenom)
pReported = p;
reportedVals = {F p};
end



function tf = redundantTerm(Q,y,inclRows)
    if nargin>=3
        y = y(inclRows,:);
    end
    yfit = Q*(Q'*y);
    res = y - yfit;
    ratio = sum(abs(res(:)).^2) / sum(abs(y(:).^2));
    tf = ratio < (eps^(3/4)) * size(Q,2) * sqrt(size(Q,1));
end

                
function predVars = getValidPredictors(X,varNames,responseVar)
% Get dataset/table variables that could be valid predictors

% Restrict predictors to specified var names
if isa(X,'dataset')
    X = dataset2tabel(X);
end
dsVarNames = X.Properties.VariableNames;
if isempty(varNames)
    predictors = 1:size(X,2);
else
    predictors = find(ismember(dsVarNames,varNames));
end

% If a response name or number is given, do not include that
if ~isempty(responseVar) && (~isnumeric(responseVar) || isscalar(responseVar))
    if islogical(responseVar)
        response = find(responseVar);
    elseif ~isnumeric(responseVar)
        response = find(strcmp(responseVar,dsVarNames));
    else
        response = responseVar;
    end
    if isscalar(response)
        predictors = predictors(predictors~=response);
    end
end

% Find the subset of variables that can be used as predictors
okResponse = false(length(predictors),1);
for j=1:length(predictors)
    v = X.(dsVarNames{predictors(j)});
    if ischar(v) && ismatrix(v)
        % okay
    elseif ~isvector(v)
        predictors(j) = 0;
    elseif iscell(v)
        if ~iscellstr(v)
            predictors(j) = 0;
        end
    elseif ~isnumeric(v) && ~islogical(v) && ~isa(v,'categorical')
        predictors(j) = 0;
    elseif isnumeric(v) || islogical(v)
        okResponse(j) = true;
    end
end

% If no response was given, find the last numbered variable that could be a
% response and remove it
if isempty(responseVar) && any(okResponse)
    predictors(find(okResponse,1,'last')) = 0;
end

% Return predictor names
predVars = dsVarNames(predictors(predictors~=0));
end

function c = countlevels(x)
if ischar(x)
    c = size(x,1);
else
    c = numel(x);
end
end
