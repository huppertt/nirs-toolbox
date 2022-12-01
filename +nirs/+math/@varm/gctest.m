function [h,Summary] = gctest(Mdl,varargin)
%GCTEST Granger causality and block exogeneity tests
%
% Syntax:
%
%   [h,Summary] = gctest(Mdl)
%   [h,Summary] = gctest(Mdl,name,value,...)
%
% Description:
%
% Granger causality and block exogeneity tests assess whether past values
% of the variable Y1 have an impact on the predictive distribution of the
% variable Y2, where Y1 and Y2 are component series in a VAR model.
%
% Supported test types:
% 1. gctest(Mdl,'Type','leave-one-out') 
% remove one variable at a time for each equation and test zero constraints
% 2. gctest(Mdl,'Type','exclude-all')
% remove all but self-lags and test Granger causality
% 3. gctest(Mdl,'Type','block-wise','Cause',...,'Effect',...)
% test zero constraints by the specified 'Cause' block and 'Effect' block
%
% Input Arguments:
%
%   Mdl -       Fitted VARM model.
%
% Optional Input Parameter Name/Value Pairs:
%
% 'Type'        String or character vector indicating test types: 
%               'leave-one-out', 'exclude-all' or 'block-wise'.
%               The default is 'leave-one-out'.
%
% 'Alpha'       Positive value of significance level for tests. 
%               The default is 0.05. 
%
% 'Test'        String or character vector indicating test statistics. 
%               Values are 'chi-square' or 'f'. The default is 'chi-square'.
%
% 'Cause'       Numeric indices, string vector or cell array of character
%               vectors indicating a subset of series for evaluating
%               one-step 'Cause'. The default is all variables.
%
% 'Effect'      Numeric indices, string vector or cell array of character 
%               vectors indicating a subset of series for evaluating 
%               one-step 'Effect'. The default is all variables.
%
% 'Display'     Logical indicator for summary table display.
%               The default is true.
%
%
% Output Arguments:
%
% h             Logical vector indicating rejection of H0. The length of
%               the vector depends on how many tests performed. Values of h
%               equal to 1 indicate rejection of the null hypothesis in
%               favor of Granger causality and endogeneity. Values of h
%               equal to 0 indicate failure to reject the null hypothesis
%               of one-step non-causality.
%
% Summary       Table that summarizes test results. Column names are: 
%               H0: null hypothesis on Granger causality/block exogeneity 
%               Decision: reject or cannot reject H0
%               Distribution: distribution of test statistics under H0
%               Statistic: value of test statistics
%               PValue: p-value of test statistics
%               CriticalValue: critical value at alpha-significant level
%
% Notes:
%
%   o Name-Value pairs 'Cause' and 'Effect' are mainly used by 'block-wise'
%     test type. In general, they are not applicable to 'leave-one-out'
%     or 'exclude-all' tests. However, if users want to run some unusual
%     tests such as constraints on self-lags, consider this:
%     gctest(Mdl,'Type','leave-one-out','Cause',1,'Effect',1); 
%     The first column of Summary (i.e., Summary.H0) shows what tests are
%     performed.
%
% References:
%
%   [1] Granger, C. W. J. Investigating Causal Relations by Econometric
%       Models and Cross-spectral Methods, Econometrica, 1969, 37, 424-459.
%
%   [2] Hamilton, J. D. Time Series Analysis. Princeton, NJ: Princeton
%       University Press, 1994.
%
%   [3] Lutkepohl, H. New Introduction to Multiple Time Series Analysis.
%       Springer-Verlag, 2007.
%
%   [4] Toda, H. Y. and T. Yamamoto. Statistical inferences in vector 
%       autoregressions with possibly integrated processes. Journal of 
%       Econometrics, 1995, 66, 225-250.
%
%   [5] Dolado, J. J. and H. Lutkepohl. Making Wald Tests Work for  
%       Cointegrated VAR Systems. Econometric Reviews, 1996, 15, 369-386.

% Copyright 2019 The MathWorks, Inc.

callerName = 'gctest';
parseObj = inputParser;
addParameter(parseObj,'Type','leave-one-out',@(x)validateattributes(x,{'char','string'},{'scalartext'},callerName));
addParameter(parseObj,'Alpha',0.05,@(x)validateattributes(x,{'numeric'},{'scalar','positive','<',1},callerName));
addParameter(parseObj,'Test','chi-square',@(x)validateattributes(x,{'char','string'},{'scalartext'},callerName));
addParameter(parseObj,'Cause',[],@(x)validateattributes(x,{'numeric','cell','char','string'},{},callerName));
addParameter(parseObj,'Effect',[],@(x)validateattributes(x,{'numeric','cell','char','string'},{},callerName));
addParameter(parseObj,'Display',true,@(x)validateattributes(x,{'numeric','logical'},{'scalar','binary'},callerName));
parse(parseObj,varargin{:});
typeStr = parseObj.Results.Type;
alpha = parseObj.Results.Alpha;
chi2Str = parseObj.Results.Test;
CauseStr = parseObj.Results.Cause;
EffectStr = parseObj.Results.Effect;
dispFlag = parseObj.Results.Display;

% Handle string inputs
if strncmpi(typeStr,'leave-one-out',1)
    typeFlag = 1;
elseif strncmpi(typeStr,'exclude-all',1)
    typeFlag = 2;
elseif strncmpi(typeStr,'block-wise',1)
    typeFlag = 3;    
else
    error(message('econ:gctest:NotSupportedType'))    
end

if strncmpi(chi2Str,'chi-square',1)
    chi2Flag = true;
elseif strncmpi(chi2Str,'f',1)
    chi2Flag = false;
else
    error(message('econ:gctest:Chi2F'))    
end

if isempty(Mdl.FitInformation)
    error(message('econ:gctest:UnfitModel'))
end

% VAR dimension
nobs = Mdl.FitInformation.SampleSize;
nvar = Mdl.NumSeries;
nlag = Mdl.P;
nX = size(Mdl.Beta,2);
numPredictors = nvar*nlag + any(Mdl.Constant~=0) + any(Mdl.Trend~=0) + nX;

% Coefficient and covariance matrix of VAR
CoeffAR = [Mdl.AR{:}];
finiteSampleAdjust = nobs / (nobs - numPredictors);
CovAR = Mdl.FitInformation.Sigma.AR .* finiteSampleAdjust;
SeriesNames = Mdl.SeriesNames;

% Enlarge the covariance matrix under coefficient constraints
% Case 1: Part of a lagged matrix is constrained, numel(CoeffAR) == size(CovAR,1)
% Case 2: All elements of a lagged matrix are constrained, numel(CoeffAR) > size(CovAR,1)
if numel(CoeffAR) > size(CovAR,1)
    isEstimated = [Mdl.FitInformation.IsEstimated.AR{:}];
    CovARSmall = CovAR;
    CovAR = zeros(numel(CoeffAR),'like',CovAR);
    CovAR(isEstimated(:),isEstimated(:)) = CovARSmall;
end

% Subset of causes and effects
CompleteFlag = isempty(CauseStr) && isempty(EffectStr);
if typeFlag == 2 && ~isempty(CauseStr)
    warning(message('econ:gctest:CauseExcludeAll'))
end
if isempty(CauseStr) && ~all(strcmp(CauseStr,''))
    Cause = 1:nvar; 
elseif isnumeric(CauseStr)
    Cause = CauseStr(:)';
    validateattributes(Cause,{'numeric'},{'positive','integer','<=',nvar},'','Cause');
else    
    Cause = find(ismember(lower(SeriesNames),lower(CauseStr)));
    Cause = Cause(:)';
    if numel(Cause) == 0
        error(message('econ:gctest:EmptyCause'))
    end
end

if isempty(EffectStr) && ~all(strcmp(EffectStr,''))
    Effect = 1:nvar;
elseif isnumeric(EffectStr)
    Effect = EffectStr(:)';
    validateattributes(Cause,{'numeric'},{'positive','integer','<=',nvar},'','Effect');
else    
    Effect = find(ismember(lower(SeriesNames),lower(EffectStr)));
    Effect = Effect(:)';
    if numel(Effect) == 0
        error(message('econ:gctest:EmptyEffect'))        
    end
end

% Summary table with columns: 
% H0, Decision, Distribution, Statistics, pValue, cValue
Summary = cell(0,6);

% Hypothesis testing of Granger causality or block exogeneity
switch typeFlag
    
    case 1
        
        %--------------
        % Leave-one-out
        %--------------
        
        % Test each variable in each equation of the VAR model
        df1 = nlag;
        df2 = nobs - numPredictors;
        for m = Effect
            for n = Cause
                
                % It is less meaningful to exclude self-lags
                % Skip such tests unless a user specifies Cause and Effect
                if m == n && CompleteFlag
                    continue
                end
                
                % Null hypothesis
                H0 = sprintf("Exclude lagged %s in %s equation",SeriesNames(n),SeriesNames(m));
                
                % Extract sub-vector and sub-matrix
                ind = m+(n-1)*nvar : nvar*nvar : nlag*nvar*nvar;
                CovARcut = CovAR(ind,ind);
                CoeffARcut = CoeffAR(ind');                
                                
                % Test statistics under the null hypothesis
                [distribution,statistics,pValue,cValue] = getTestStatistics(CovARcut,CoeffARcut,chi2Flag,alpha,df1,df2);
                                
                % Add statistics to summary table
                Summary = addSummary(Summary,statistics,pValue,cValue,alpha,distribution,H0);                
            end
        end
        
    case 2
        
        %--------------
        % Exclude-all
        %--------------
                        
        % Test for each equation in the VAR model
        df1 = nlag * (nvar-1);
        df2 = nobs - numPredictors;
        for m = Effect
            
            % Null hypothesis
            H0 = sprintf("Exclude all but lagged %s in %s equation",SeriesNames(m),SeriesNames(m));
            
            % Extract sub-vector and sub-matrix            
            ind = m : nvar: nlag*nvar*nvar;
            ind(m:nvar:end) = [];
            CovARcut = CovAR(ind,ind);
            CoeffARcut = CoeffAR(ind');  
                        
            % Test statistics under the null hypothesis
            [distribution,statistics,pValue,cValue] = getTestStatistics(CovARcut,CoeffARcut,chi2Flag,alpha,df1,df2);
            
            % Add statistics to summary table
            Summary = addSummary(Summary,statistics,pValue,cValue,alpha,distribution,H0);
        end
        
    case 3
        
        %--------------
        % Block-wise
        %--------------
        
        % Extract sub-vector and sub-matrix
        IndMat = reshape(1:nvar*nvar*nlag,[nvar,nvar,nlag]);
        ind = zeros(0,1);
        for m = 1:nlag
            submatrix = IndMat(Effect,Cause,m);
            ind = [ind;submatrix(:)]; %#ok<AGROW>
        end        
        CovARcut = CovAR(ind,ind);
        CoeffARcut = CoeffAR(ind);
        
        % Test statistics under the null hypothesis
        % For F(df1,df2) test, assume that df2 = nobs - numPredictors
        % It is also reasonable to use df2 = (nobs - numPredictors) * nvar
        df1 = numel(CoeffARcut);
        df2 = nobs - numPredictors;
        [distribution,statistics,pValue,cValue] = getTestStatistics(CovARcut,CoeffARcut,chi2Flag,alpha,df1,df2);
        
        % Add statistics to summary table
        H0 = sprintf("Exclude lagged %s in %s equations",strjoin(SeriesNames(Cause),','),strjoin(SeriesNames(Effect),','));
        Summary = addSummary(Summary,statistics,pValue,cValue,alpha,distribution,H0);
    
end

% Convert cell to table
Summary = cell2table(Summary,'VariableNames',{'H0','Decision','Distribution','Statistic','PValue','CriticalValue'});
h = strcmp(Summary.Decision,"Reject H0");
if dispFlag
    disp(Summary)
end

end


%-------------------------------------------------------------------------
% Add a new entry (row) to the summary table
function Summary = addSummary(Summary,statistics,pValue,cValue,alpha,distribution,H0)
if pValue <= alpha
    decision = "Reject H0";
else
    decision = "Cannot reject H0";
end
count = size(Summary,1);
count = count + 1;
Summary(count,:) = {H0,decision,distribution,statistics,pValue,cValue};
end


%-------------------------------------------------------------------------
% Wald chi2 or F test statistics
% A fitted VARM may be subject to zero constraints in AR coefficients
function [distribution,statistics,pValue,cValue] = getTestStatistics(CovARcut,CoeffARcut,chi2Flag,alpha,df1,df2)
mask = diag(CovARcut) == 0;
if any(mask) && all(CoeffARcut(mask)==0)
    CovARcut = CovARcut(~mask,~mask);
    CoeffARcut = CoeffARcut(~mask);
    df1 = df1 - sum(mask);
end
df1(df1==0) = 0.001;
[CovARChol,flag] = chol(CovARcut,'lower');
if flag > 0
    CovARcut = CovARcut + 1e-8 * eye(size(CovARcut),'like',CovARcut);
    CovARChol = chol(CovARcut,'lower');
end
normalRV = CovARChol \ CoeffARcut;
if chi2Flag
    distribution = sprintf("Chi2(%d)",df1);
    statistics = normalRV' * normalRV;
    pValue = chi2cdf(statistics,df1,'upper');
    cValue = chi2inv(1-alpha,df1);
else
    distribution = sprintf("F(%d,%d)",df1,df2);
    statistics = normalRV' * normalRV / df1;
    pValue = fcdf(statistics,df1,df2,'upper');
    cValue = finv(1-alpha,df1,df2);
end
end

