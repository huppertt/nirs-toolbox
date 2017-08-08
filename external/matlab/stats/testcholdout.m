function [h,p,err1,err2] = testcholdout(Yhat1,Yhat2,Y,varargin)
%TESTCHOLDOUT Compare accuracies of two sets of predicted class labels.
%   H=TESTCHOLDOUT(YHAT1,YHAT2,Y) performs a test of the null hypothesis: two
%   sets of class labels, YHAT1 and YHAT2, have equal accuracy for predicting
%   true class labels Y.
%       H = 0 => Do not reject the null hypothesis at the 5% significance level.
%       H = 1 => Reject the null hypothesis at the 5% significance level.
%
%   Pass YHAT1, YHAT2 and Y as categorical arrays, character arrays, logical
%   vectors, numeric vectors, or cell arrays of strings. If you pass these as
%   character arrays, they must have one class label per row. YHAT1, YHAT2 and Y
%   must have the same number of elements.
%
%   [H,P,E1,E2]=TESTCHOLDOUT(YHAT1,YHAT2,Y) also returns the p-value P of the
%   test, classification error E1 for YHAT1 and classification error E2 for
%   YHAT2. If you pass the 'Cost' parameter, E1 and E2 are misclassification
%   costs.
%
%   [...]=TESTCHOLDOUT(YHAT1,YHAT2,Y,'PARAM1',val1,'PARAM2',val2,...) specifies
%   one or more of the following name/value pairs:
%       'Alpha'        - Confidence level, a positive scalar. Default: 0.05
%       'Alternative'  - String indicating the alternative hypothesis, one of:
%                          * 'unequal' - TESTCHOLDOUT tests H0: "YHAT1 and YHAT2
%                                        have equal accuracy" against H1: "YHAT1
%                                        and YHAT2 have unequal accuracy".
%                          * 'less'    - TESTCHOLDOUT tests H0: "YHAT1 is at
%                                        least as accurate as YHAT2" against H1:
%                                        "YHAT1 is less accurate than YHAT2".
%                          * 'greater' - TESTCHOLDOUT tests H0: "YHAT1 is at
%                                        most as accurate as YHAT2" against H1:
%                                        "YHAT1 is more accurate than YHAT2".
%                          Default: 'unequal'
%       'ClassNames'   - Array of class names. Use the data type that exists in
%                        Y. You can use this argument to order the classes or
%                        select a subset of classes. Default: The class names
%                        that exist in Y.
%       'Cost'         - Square matrix, where COST(I,J) is the cost of
%                        classifying a point into class J if its true class is
%                        I. Alternatively, COST can be a structure S having two
%                        fields: S.ClassificationCosts containing the cost
%                        matrix C, and S.ClassNames containing the class names
%                        and defining the ordering of classes used for the rows
%                        and columns of the cost matrix. For S.ClassNames use
%                        the data type that exists in Y. If you pass 'Cost' as a
%                        numeric matrix, the order of rows and columns matches
%                        the order defined by 'ClassNames'. If you pass 'Cost',
%                        TESTCHOLDOUT can perform an asymptotic two-sided test
%                        only, that is, the 'Alternative' and 'Test' parameters
%                        must be set to 'unequal' and 'asymptotic',
%                        respectively. Default: []
%       'CostTest'     - String indicating the type of the cost-sensitive test,
%                        one of: 'likelihood' or 'chisquare'. If 'likelihood',
%                        TESTCHOLDOUT uses a likelihood ratio test. If
%                        'chisquare', TESTCHOLDOUT uses a chi-square test. This
%                        parameter is ignored unless you pass a cost matrix
%                        using the 'Cost' parameter. The chi-square test
%                        requires an Optimization Toolbox license. Type 'doc
%                        testcholdout' for info about cost-sensitive tests.
%       'Test'         - String, one of: 'asymptotic', 'exact', or 'midp'. Let
%                        N01 be the number of observations misclassified by
%                        YHAT1 and correctly classified by YHAT2. Let N10 be the
%                        number of observations correctly classified by YHAT1
%                        and misclassified by YHAT2. If you set 'Test' to
%                          * 'asymptotic' - TESTCHOLDOUT can perform several
%                                           tests:
%                               o If you do not pass 'Cost', TESTCHOLDOUT
%                                 performs the asymptotic McNemar test assuming
%                                 (N01-N10)/sqrt(N01+N10) has a normal
%                                 distribution with zero mean and unit variance.
%                               o If you pass 'Cost', TESTCHOLDOUT by default
%                                 performs a likelihood ratio test assuming
%                                 twice the log of the likelihood ratio has a
%                                 chi-square distribution with one degree of
%                                 freedom. If you set 'CostTest' to 'chisquare',
%                                 TESTCHOLDOUT performs a chi-square test
%                                 assuming the test statistic has a chi-square
%                                 distribution with one degree of freedom.
%                          * 'exact'      - TESTCHOLDOUT performs the exact
%                                           conditional McNemar test assuming
%                                           that N01 has a binomial distribution
%                                           with N01+N10 trials and binomial
%                                           parameter 0.5.
%                          * 'midp'       - TESTCHOLDOUT performs the mid-p
%                                           value McNemar test. This test uses
%                                           the same distribution assumptions as
%                                           the exact conditional test does. To
%                                           compute the p-value, this test
%                                           corrects the binomial CDF
%                                           P(X<=x;N01+N10,0.5) by subtracting
%                                           half of the binomial probability
%                                           mass function at x,
%                                           0.5*P(X=x;N01+N10,0.5).
%                          Default: 'midp' if 'Cost' is not passed and
%                                   'asymptotic' otherwise
%
% Example: Compare SVM and bagged trees on ionosphere data using a held-out set. 
%   load ionosphere;
%   rng(151515); % set the RNG seed for reproducibility
%   cvp = cvpartition(Y,'holdout',0.5);
%   itrain = training(cvp);
%   itest = test(cvp);
%   svm = fitcsvm(X(itrain,:),Y(itrain),'Standardize',true,...
%       'KernelFunction','RBF','KernelScale','auto');
%   YhatSVM = predict(svm,X(itest,:));
%   bag = fitensemble(X(itrain,:),Y(itrain),'Bag',100,'Tree','type','classification');
%   YhatBag = predict(bag,X(itest,:));
%   [h,p] = testcholdout(YhatSVM,YhatBag,Y(itest))
%
%   See also testckfold.

%   Copyright 2014 The MathWorks, Inc.

Y = classreg.learning.internal.ClassLabel(Y);
Y1 = classreg.learning.internal.ClassLabel(Yhat1);
Y2 = classreg.learning.internal.ClassLabel(Yhat2);
nonzeroClassNames = levels(Y); 

% Check the sizes.
N1 = numel(Y1);
if numel(Y)~=N1
    error(message('stats:testcholdout:ClassLabelSizeMismatch','Yhat1'));
end

N2 = numel(Y2);
if N1~=N2
    error(message('stats:testcholdout:ClassLabelSizeMismatch','Yhat2'));
end

% Decode input args
args = {'alpha' 'alternative' 'test' 'classnames' 'cost'   'costtest'};
defs = {   0.05     'unequal'     ''           ''     [] 'likelihood'};

[alpha,alternative,mode,userClassNames,cost,costtest] = ...
    internal.stats.parseArgs(args,defs,varargin{:});

% Check alpha
if ~isscalar(alpha) || ~isfloat(alpha) || ~isreal(alpha) || isnan(alpha) ...
        || alpha<=0 || alpha>=1
    error(message('stats:testcholdout:BadAlpha'));
end

% Is this a cost-sensitive test?
fitcost = false;
if ~isempty(cost)
    fitcost = true;
end

% Determine the alternative hypothesis.
alternative = validatestring(alternative,{'unequal' 'less' 'greater'},...
    'testcholdout','Alternative');

% Determine the test type. Make sure the test type is compatible with
% cost-sensitive analysis.
if isempty(mode) 
    if fitcost
        mode = 'asymptotic';
    else
        mode = 'midp';
    end
end

mode = validatestring(mode,{'asymptotic' 'exact' 'midp'},...
    'testcholdout','Test');

if ~(strcmp(mode,'asymptotic') && strcmp(alternative,'unequal')) && fitcost
    error(message('stats:testcholdout:BadCostAlternativeTestCombo'));
end

% Determine the type of the cost-sensitive test
costtest = validatestring(costtest,{'likelihood' 'chisquare'},...
    'testcholdout','CostTest');

% Remove rows for which true class values are missing
t = ismissing(Y);
if all(t)
    error(message('stats:testcholdout:AllObservationsHaveMissingTrueLabels'));
end
if any(t)
    Y(t)  = [];
    Y1(t) = [];
    Y2(t) = [];
end

% Process class names.
if isempty(userClassNames)
    % If the user has not passed any class names, use those found in the array
    % of true class labels.
    userClassNames = nonzeroClassNames;
else
    userClassNames = classreg.learning.internal.ClassLabel(userClassNames);
    
    % If none of the class names passed by the user is found in the existing
    % class names, error.
    missingC = ~ismember(userClassNames,nonzeroClassNames);
    if all(missingC)
        error(message('stats:classreg:learning:classif:FullClassificationModel:prepareData:ClassNamesNotFound'));
    end
    
    % If the user passed a subset of classes found in the data, remove labels
    % for classes not included in that subset.
    missingC = ~ismember(nonzeroClassNames,userClassNames);
    if any(missingC)
        unmatchedY = ismember(Y,nonzeroClassNames(missingC));
        Y(unmatchedY)  = [];
        Y1(unmatchedY) = [];
        Y2(unmatchedY) = [];
        nonzeroClassNames(missingC) = [];
    end
end

% Record the number of observations and useful classes.
N = numel(Y);
K = numel(nonzeroClassNames);

if fitcost % cost-sensitive analysis
    % Remove entries for useless classes and check the validity of the cost
    % matrix.
    cost = classreg.learning.classif.FullClassificationModel.processCost(...
        cost,ones(1,K),userClassNames,nonzeroClassNames);
    
    % Logical matrix of size N-by-K for the true class membership.
    C  = classreg.learning.internal.classCount(nonzeroClassNames,Y);
    
    % Prepare C1 and C2, similar to C. We cannot use the classCount function
    % here because it errors when the 2nd input has levels not present in the
    % 1st input. Y1 and Y2 can contain classes not found in Y.
    C1 = false(N,K);
    for k=1:K
        C1(:,k) = ismember(Y1,nonzeroClassNames(k));
    end
    
    C2 = false(N,K);
    for k=1:K
        C2(:,k) = ismember(Y2,nonzeroClassNames(k));
    end
    
    % Compute misclassification costs.
    err1 = mean(sum((C*cost).*C1,2));
    err2 = mean(sum((C*cost).*C2,2));
    
    % Form a K-by-K-by-K matrix of observed cell counts.
    Mobs = zeros(K,K,K);
    for i=1:K
        for j=1:K
            for k=1:K
                Mobs(i,j,k) = sum( Y1==nonzeroClassNames(i) ...
                    & Y2==nonzeroClassNames(j) & Y==nonzeroClassNames(k) );
            end
        end
    end
    
    % Fit cell probabilities
    [~,~,p] = fitCellCounts(Mobs,cost,costtest);
    
else % regular (cost-insensitive) analysis
    
    % Count correct predictions
    good1 = Y1==Y;
    good2 = Y2==Y;
    
    % Compute misclassification errors
    err1 = mean(~good1);
    err2 = mean(~good2);
    
    % Compute contingency table cells
    n01 = sum(~good1 & good2); % Missed by C1 and found by C2
    n10 = sum(good1 & ~good2); % Found by C1 and missed by C2

    % Edge cases
    
    if n01==0 && n10==0
        h = false;
        p = 1;
        return;
    end
    
    if n01==n10 && strcmp(alternative,'unequal')
        h = false;
        p = 1;
        return;
    end

    % Perform the test
    switch mode
        case 'asymptotic'
            switch alternative
                case 'unequal'
                    p = chi2cdf((n01-n10)^2/(n01+n10),1,'upper');
                case 'less'
                    p = normcdf((n10-n01)/sqrt(n01+n10)); % large n01-n10
                case 'greater'
                    p = normcdf((n01-n10)/sqrt(n01+n10)); % small n01-n10
            end
            
        case 'exact'
            switch alternative
                case 'unequal'
                    p = 2*binocdf(min(n01,n10),n01+n10,0.5);
                case 'less'
                    p = binocdf(n10,n01+n10,0.5); % large n01 and small n10
                case 'greater'
                    p = binocdf(n01,n01+n10,0.5); % small n01 and large n10
            end
            
        case 'midp'
            switch alternative
                case 'unequal'
                    p = 2*binocdf(min(n01,n10)-1,n01+n10,0.5) + ...
                        binopdf(min(n01,n10),n01+n10,0.5);
                case 'less'
                    % large n01 and small n10
                    p = binocdf(n10-1,n01+n10,0.5) + 0.5*binopdf(n10,n01+n10,0.5); 
                case 'greater'
                    % small n01 and large n10
                    p = binocdf(n01-1,n01+n10,0.5) + 0.5*binopdf(n01,n01+n10,0.5); 
            end
            
    end
    
    
end

% Accept = 0, reject = 1
h = p<alpha;

end


function [Mhat,chisqval,pval] = fitCellCounts(Mobs,C,costtest)
%fitCellCounts Find cell counts expected under the null hypothesis of equal cost.
%   [Mhat,chisqval,pval]=fitCellCounts(Mobs,C,costtest) returns expected cell
%   counts Mhat, the value of the test statistic chisqval and the associated
%   p-value pval for observed cell counts Mobs and cost matrix C. Pass Mobs as a
%   K-by-K-by-K array of non-negative integers. Pass C as a K-by-K matrix of
%   non-negative values with zeros on the main diagonal. Pass costtest as a
%   string, one of: 'likelihood' or 'chisquare'. The test statistic has an
%   asymptotic chi-square distribution with 1 d.o.f.

% Get the number of classes
K = size(C,1);

% Check the cell count matrix
if ~isfloat(Mobs) || (ndims(Mobs)~=3 && ~isscalar(Mobs)) || any(size(Mobs)~=K) ...
        || any(Mobs(:)<0) || any(round(Mobs(:))~=Mobs(:))
    error(message('stats:testcholdout:BadMatrixOfCellCounts',K,K,K));
end

% Check the matrix of misclassification costs
if ~isfloat(C) || ~ismatrix(C) || any(size(C)~=K) ...
        || any(C(1:K+1:end)~=0) || any(C(:)<0)
    error(message('stats:testcholdout:BadCostMatrix'));
end

% Edge case: misclassification costs nothing
if all(C(:)==0)
    Mhat = Mobs + 1;
    chisqval = 0;
    pval = 1;
    return;
end

% Flatten out the cost matrix in two arrays representing the cost of the 1st and
% 2nd classifiers. We do this to be able to impose the linear constraint.
cost1 = zeros(K^3,1);
cost2 = zeros(K^3,1);

for k=1:K
    a1 = repmat(C(k,:)',1,K);
    a2 = repmat(C(k,:) ,K,1);
    idx = (k-1)*K^2+1 : k*K^2;
    cost1(idx) = a1(:);
    cost2(idx) = a2(:);
end

% Vector of cost differences, c_{ki}-c_{kj}. Use it as Aeq to impose the
% equality constraint for quadprog.
dcost = cost1 - cost2;

% Flatten out the array of cell counts
Mobs = Mobs(:);

% Laplace correction
Mobs = Mobs+1; 

%
% Solve for optimal cell counts: Mhat(i,j,k) = sum(Mobs(:))*p(i,j,k) for cell
% probability p(i,j,k)
%

switch costtest 
    case 'likelihood' % likelihood ratio test
        
        N = sum(Mobs);
        
        % Set limits on lambda. Shrink the theoretical limits a bit to prevent
        % fzero from overflowing at the endpoints.
        maxdcost = max(abs(dcost));
        maxlambda = N/maxdcost; % theoretical endpoints
        maxlambda = maxlambda*(1 - eps(class(maxlambda)));
        lims = [-maxlambda maxlambda];
        
        % Find the optimal lambda under H0
        f = @(lambda) Mobs'*(dcost./(N + lambda*dcost));
        try
            lambda = fzero(f,lims);
        catch me
            error(message('stats:testcholdout:FzeroFails',me.message));
        end
        
        % Estimated cell counts under H0. (Estimated cell counts under H1 are Mobs.)
        lamdcost = lambda*dcost/N;
        Mhat = Mobs./(1+lamdcost);
        
        % Log likelihood ratio
        chisqval = 2*Mobs'*log1p(lamdcost);
        
    case 'chisquare' % Use quadprog
        
        % Check the Optim license
        if isempty(ver('Optim'))
             error(message('stats:testcholdout:NeedOptimForChisquareTest'));
        end
        
        % Equality constraint: sum p_{ijk}*(c_{ki}-c_{kj} = 0
        beq = 0;
        
        % Prepare an indicator vector to exclude cells {ijk} with i==j from fitting.
        indicator = ones(K^3,1);
        
        for k=1:K
            ibegin = (k-1)*K^2+1;
            iend   = k*K^2;
            indicator(ibegin:K+1:iend) = 0;
        end
        
        H = diag(sparse(2*indicator./Mobs)); % Hessian for quadprog
        
        f = -2*indicator; % linear part for quadprog
        
        %opts = optimoptions(@quadprog,'Display','iter-detailed','Diagnostics','on');
        opts = optimoptions(@quadprog,'Display','none');
        
        % Solve
        [Mhat,~,exitflag] = quadprog(H,f,[],[],dcost',beq,[],[],[],opts);
        
        if exitflag~=1
            warning(message('stats:testcholdout:BadExitFlagFromQuadprog',exitflag));
        end
        
        % Compute the chi-square value
        chisqval = sum( indicator .* (Mhat-Mobs).^2 ./ Mobs );
end
    
% Reshape the solution vector back into a matrix of cell counts
Mhat = reshape(Mhat,K,K,K);

% Get the p-value
pval = chi2cdf(chisqval,1,'upper');

end

