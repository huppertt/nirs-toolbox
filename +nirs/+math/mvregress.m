function [Param,Covar,Resid,VarParam,Objective] = mvregress(Design, Data, varargin)
%MVREGRESS Multivariate regression with missing data.
%   [BETA,SIGMA,RESID]=MVREGRESS(X,Y) performs multivariate regression of
%   the multivariate observations in the N-by-D matrix Y on the predictor
%   variables in X, and returns a P-by-1 column vector BETA of coefficient
%   estimates, a D-by-D matrix SIGMA of the estimated covariance of Y, and an
%   N-by-D matrix RESID of residuals.  NaN values in X or Y are taken to be
%   missing.  Observations with missing values in X are ignored.  Missing
%   values in Y are handled according to the value of the 'algorithm'
%   parameter described below.
%
%   For any value of D>1 with X being an N-by-P matrix, it returns a 
%   P-by-D matrix BETA of coefficient estimates.
%
%   Y is an N-by-D matrix of D-dimensional multivariate observations.  X may
%   be either a matrix or a cell array.  If D=1, X may be an N-by-P design
%   matrix of predictor variables.  For any value of D, X may also be a cell
%   array of length N, each cell containing a D-by-P design matrix for one
%   multivariate observation.  If all observations have the same D-by-P
%   design matrix, X may be a single cell. For any value of D, X may also
%   be an N-by-P matrix. 
%
%   [BETA,SIGMA,RESID,VARPARAM]=MVREGRESS(...) also returns an estimated
%   covariance matrix of the estimates.  By default, or if the 'varformat'
%   parameter is 'beta' (see below), VARPARAM is the estimated covariance
%   matrix of the coefficient estimates BETA.  If the 'varformat' parameter
%   is 'full', VARPARAM is the estimated covariance matrix of the combined
%   BETA and SIGMA estimates.
%
%   [BETA,SIGMA,RESID,VARPARAM,OBJECTIVE]=MVREGRESS(...) also returns the
%   value of the objective function, or log likelihood, after the last
%   iteration.
%
%   [...]= MVREGRESS(X,Y,'PARAM1',VALUE1,'PARAM2',VALUE2,...) specifies
%   additional parameter name/value pairs chosen from the following:
%     'algorithm'  Either 'ecm' to compute the maximum likelihood estimates
%                  via the ECM algorithm, 'cwls' to perform least squares
%                  (optionally conditionally weighted by an input covariance
%                  matrix), or 'mvn' to omit observations with missing data
%                  and compute the ordinary multivariate normal estimates.
%                  Default is 'mvn' for complete data, 'ecm' for missing data
%                  when the sample size is sufficient to estimate all
%                  parameters, and 'cwls' otherwise.
%     'covtype'    Either 'full' (default) to allow a full covariance matrix,
%                  or 'diagonal' to restrict it to be a diagonal matrix.
%     'maxiter'    Maximum number of iterations (default 100).
%     'tolbeta'    Convergence tolerance for BETA (default sqrt(eps)).
%                  Iterations continue until the TOLBETA and TOLOBJ conditions
%                  are met.  The test for convergence at iteration k is
%                     ||BETA(k)-BETA(k-1)|| < sqrt(P)*TOLBETA * (1+||BETA(k)||)
%                  where ||v|| represents the norm of the vector v.
%     'tolobj'     Convergence tolerance for changes in the objective function
%                  (default eps^(3/4)).  The test is
%                     |Obj(k)-Obj(k-1)| < TolObj * (1 + |Obj(k)|)
%                  If both TOLOBJ and TOLBETA are 0, the function performs
%                  MAXITER iterations with no convergence test.
%     'beta0'      A vector of P elements to be used as the initial estimate
%                  for BETA.  Default is a zero vector.  For any value of D>1 
%                  with X being an N-by-P matrix, beta0 should be a P-by-D 
%                  matrix, default is a zero N-by-P matrix. Not used for 
%                  the 'mvn' algorithm.
%     'covar0'     A D-by-D matrix to be used as the initial estimate for
%                  SIGMA.  Default is the identity matrix.  For the 'cwls'
%                  algorithm, this matrix is usually a diagonal matrix and it
%                  is not changed during the iterations, so the input value
%                  is used as the weighting matrix at each iteration.
%     'outputfcn'  An output function.
%     'varformat'  Either 'beta' to compute VARPARAM for BETA only (default),
%                  or 'full' to compute VARPARAM for both BETA and SIGMA.
%     'vartype'    Either 'hessian' to compute VARPARAM using the Hessian or
%                  observed information (default), or 'fisher' to compute the
%                  complete-data Fisher or expected information.  The 'hessian'
%                  method takes into account the increased uncertainties due
%                  to missing data, while the 'fisher' method does not.
%
%   The RESID values corresponding to missing values in Y are the differences
%   between the conditionally-imputed values for Y and the fitted values.  The
%   SIGMA estimate is not the sample covariance matrix of the RESID matrix.
%
%   The output function is called with three arguments:
%      1.  Vector of current parameter estimates
%      2.  A structure with fields 'Covar' for the current value of the
%          covariance matrix, 'iteration' for the current iteration number,
%          and 'fval' for the current value of the objective function.
%      3.  A text string that is 'init' when called during initialization,
%          'iter' when called after an iteration, and 'done' when called
%          after completion.
%   The function should return TRUE if the iterations should stop, or FALSE
%   if they should continue.
%
%   Example: Predict regional flu estimates based on Google queries using
%            the national CDC estimates as a predictor
%      load flu
%      y = double(flu(:,2:end-1));  % response = regional queries
%      x = flu.WtdILI;              % predictor = national CDC estimates
%      [nobs,nregions] = size(y);
%
%      % Create and fit model with separate intercepts but common slope
%      X = cell(nobs,1);
%      for j=1:nobs
%         X{j} = [eye(nregions), repmat(x(j),nregions,1)];
%      end
%      [b,sig,resid,vars,loglik] = mvregress(X,y);
%
%      % Plot raw data with fitted lines
%      B = [b(1:nregions)';repmat(b(end),1,nregions)]
%      subplot(3,1,1);
%      xx = linspace(.5,3.5)';
%      h = plot(x,y,'x', xx, [ones(size(xx)),xx]*B,'-');
%      for j=1:nregions; set(h(nregions+j),'color',get(h(j),'color')); end
%      regions = flu.Properties.VarNames;
%      legend(regions{2:end-1},'location','NW')
%
%      % Create and fit model with separate intercepts and slopes
%      for j=1:nobs
%         X{j} = [eye(nregions), x(j)*eye(nregions)];
%      end
%      [b,sig,resid,vars,loglik2] = mvregress(X,y);
%
%      % Plot raw data with fitted lines
%      B = [b(1:nregions)';b(nregions+1:end)']
%      subplot(3,1,2);
%      h = plot(x,y,'x', xx, [ones(size(xx)),xx]*B,'-');
%      for j=1:nregions; set(h(nregions+j),'color',get(h(j),'color')); end
%
%      % Likelihood ratio test for significant difference
%      chisq = 2*(loglik2-loglik)
%      p = 1-chi2cdf(chisq, nregions-1)
%
%      % Create and fit model with separate intercepts and slopes in matrix form
%      X = [ones(size(x)),x];
%      [b,sig,resid,vars,loglik2] = mvregress(X,y);
%
%      % Plot raw data with fitted lines
%      B = b
%      subplot(3,1,3);
%      h = plot(x,y,'x', xx, [ones(size(xx)),xx]*B,'-');
%      for j=1:nregions; set(h(nregions+j),'color',get(h(j),'color')); end
%
%   See also REGSTATS, MANOVA1, MVREGRESSLIKE.

% References:
%    [1] Roderick J. A. Little and Donald B. Rubin, Statistical Analysis with
%        Missing Data, 2nd ed., John Wiley & Sons, Inc., 2002.
%    [2] Xiao-Li Meng and Donald B. Rubin, "Maximum Likelihood Estimation via
%        the ECM Algorithm," Biometrika, Vol. 80, No. 2, 1993, pp. 267-278.
%    [3] Joe Sexton and Anders Rygh Swensen, "ECM Algorithms that Converge at
%        the Rate of EM," Biometrika, Vol. 87, No. 3, 2000, pp. 651-662.
%    [4] A. P. Dempster, N.M. Laird, and D. B. Rubin, "Maximum Likelihood from
%        Incomplete Data via the EM Algorithm," Journal of the Royal Statistical
%        Society, Series B, Vol. 39, No. 1, 1977, pp. 1-37.

%    Copyright 2006-2017 The MathWorks, Inc.


% Step 1 - check arguments
if nargin > 2
    [varargin{:}] = convertStringsToChars(varargin{:});
end

if nargin < 2
    error(message('stats:mvregress:MissingInputArg'));
end
if isempty(Data)
    error(message('stats:statcheckmvnr:EmptyDataArray'));
end

if isempty(Design)
    error(message('stats:statcheckmvnr:EmptyDesignArray'));
end

okargs =   {'maxiter'   'tolparam'  'tolobj'    'param0'  'covar0'  'beta0' ...
            'algorithm' 'outputfcn' 'varformat' 'vartype' 'covtype' 'tolbeta'};
defaults = {100         sqrt(eps)   eps^(3/4)   []        []        [] ...
            []          []          'beta'      'hessian' 'full'    sqrt(eps)};
[           MaxIter,     TolParam,    TolObj,      Param0,    Covar0,    Beta0, ...
            EstMethod,   OutFun,      VarFormat,   VarType,   CovType,   TolBeta] = ...
                internal.stats.parseArgs(okargs,defaults,varargin{:});

% Some options have synonyms for compatibility. Take the documented value
% if both are given.
if ~isempty(Beta0)
    Param0 = Beta0;
end
if ~isempty(TolBeta)
    TolParam = TolBeta;
end

if ~isscalar(MaxIter)
    error(message('stats:mvregress:BadMaxIter'));
elseif MaxIter < 1
    MaxIter = 1;
end

if ~isempty(EstMethod)
    if ~ischar(EstMethod) || size(EstMethod,1)~=1
        EstMethod = [];
    else
        [~,EstMethod] = internal.stats.getParamVal(EstMethod,{'cwls' 'ecm' 'mvn'},'algorithm');
    end
    if isempty(EstMethod)
        error(message('stats:mvregress:BadAlgorithm'));
    end
end
if ~(isempty(EstMethod) || ismember(EstMethod,1:3))
    error(message('stats:mvregress:BadAlgorithm'));
end

% Check variance format
okvals = {'beta' 'full'};
if ~ischar(VarFormat) || size(VarFormat,1)~=1
    VarFormat = [];
else
    [~,VarFormat] = internal.stats.getParamVal(VarFormat,okvals,'varformat');
end
if isempty(VarFormat)
    error(message('stats:mvregress:BadVarFormat'));
end
okvals{1} = 'paramonly';  % internal code for 'beta'
VarFormat = okvals{VarFormat};

% Check variance type
okvals = {'hessian' 'fisher'};
if ~ischar(VarType) || size(VarType,1)~=1
    VarType = [];
else
    [~,VarType] = internal.stats.getParamVal(VarType,okvals,'vartype');
end
if isempty(VarType)
    error(message('stats:mvregress:BadVarType'));
end
VarType = okvals{VarType};

% Check covariance type
okvals = {'full' 'diagonal'};
if ~ischar(CovType) || size(CovType,1)~=1
    CovType = [];
else
    [~,CovType] = internal.stats.getParamVal(CovType,okvals,'covtype');
end
if isempty(CovType)
    error(message('stats:mvregress:BadCovType'));
end
isdiagonal = isequal(okvals{CovType},'diagonal');

% Check inputs, ignoring NaN rows for mvn method (3)
if isempty(EstMethod) && ~any(isnan(Data(:)))
    EstMethod = 3;   % use faster method for cases where others are equivalent
end
[NumSamples, NumSeries, NumParams, Data, Design, goodrows] = ...
          statcheckmvnr(Data, Design, Param0, Covar0, isequal(EstMethod,3));

celldesign = iscell(Design);
if celldesign && (numel(Design) == 1)
    SingleDesign = true;
else
    SingleDesign = false;
end

% Step 2 - observability and ignorability tests
Count = sum(all(isnan(Data),2));
if ((NumSamples - Count) * NumSeries) <= max(NumParams, (NumSeries * (NumSeries + 1))/2)
    if ((NumSamples - Count) * NumSeries) <= NumParams
        error(message('stats:mvregress:InsufficientData'));
    elseif isempty(EstMethod) || EstMethod ~= 1  % other than cwls method
        if isempty(EstMethod)
            EstMethod = 1; % select cwls
        else
            warning(message('stats:mvregress:TryingReducedModel'));
            if EstMethod == 2  % ecm method
                EstMethod = 1; % switch to cwls
            else
                MaxIter = 1;   % or for mvn method, do just 1 iteration
            end
        end
    end
end
if isempty(EstMethod)
    EstMethod = 2;  % ecm method
end

% Step 3 - initialization
if isempty(Param0)
    Param = zeros(NumParams, 1);
else
    Param = Param0;
end

if isempty(Covar0)        % setup covariance-weighted least-squares
    CovarCWLS = eye(NumSeries, NumSeries);
    Covar = CovarCWLS;
    cwlsdiagonal = true;
else
    CovarCWLS = Covar0;
    Covar = Covar0;
    if EstMethod == 1
        cwlsdiagonal = isequal(Covar, diag(diag(Covar)));
    end
end
if isdiagonal
    Covar = diag(diag(Covar));
end

VarParam = [];

[CholCovarCWLS, CholState] = chol(CovarCWLS);
if CholState > 0
    warning(message('stats:mvregress:NonPosDefCovar'));
    CovarCWLS = eye(NumSeries, NumSeries);
    CholCovarCWLS = eye(NumSeries, NumSeries);
end

Resid = nan(NumSamples, NumSeries);

Design0 = Design;
if ~SingleDesign
    if iscell(Design)
       % Xdesign = reshape(permute(cat(3,Design{:}),[1 3 2]),[NumSeries*NumSamples,NumParams]);
        for i=1:length(Design)
            Design{i}=Design{i}';
        end

        Xdesign = cat(2,Design{:})';

    else
        % X is N-by-P matrix with D>1, result is P-by-D 
        c = cell(NumSamples,1);
        for j = 1:NumSamples
            c{j} = kron(eye(NumSeries),Design(j,:));
        end
        Design = c;
        NumParams0 = NumParams;
        [~, NumParams] = size(Design{1});
        if isempty(Param0)
            Param = zeros(NumParams,1);
        else
            Param = reshape(Param0,[NumParams 1]);
        end
        Xdesign = reshape(permute(cat(3,Design{:}),[1 3 2]),[NumSeries*NumSamples,NumParams]);
    end
else
    Xdesign = repmat(Design{1},NumSamples,1);
end

if ~isempty(OutFun)
    str.Covar = CovarCWLS;
    str.iteration = 0;
    str.fval = [];
    if OutFun(Param,str,'init')
        disp(getString(message('stats:mvregress:TerminatedByUser')));
        Resid = [];
        return
    end
end

% Step 4 - main loop
nans = isnan(Data);
partialrows = find(any(nans,2))';  % these rows already removed for mvn method
Count = sum(~all(isnan(Data),2));
seps = sqrt(eps);

% Initialize this variable so it will not appear to be a function here
Objective0 = 0;

for Iter = 1:MaxIter
    Z = Data;
    if ~isempty(partialrows)  % always empty for mvn method
        if celldesign && SingleDesign
            Mean = Design{1} * Param;
        else
            Means = reshape(Xdesign*Param, [NumSeries, NumSamples]);
        end

        % Step 5 - parameter combined E and CM step
        WarnState = warning('off','MATLAB:nearlySingularMatrix');
        for i = partialrows
            if ~SingleDesign
                Mean = Means(:,i);
            end

            [mX, mY, ~, CXY, CYY] = ecmpart(Data(i,:), Mean, Covar);

            if ~isempty(mY)
                P = isnan(Data(i,:));
                Q = ~P;
                Y = Data(i,Q)';

                Z(i,P) = mX + CXY * (CYY \ (Y - mY));
                Z(i,Q) = Y;
            end
        end
        warning(WarnState);
    end

    A = reshape(CholCovarCWLS' \ reshape(Xdesign,NumSeries,NumSamples*NumParams),[NumSeries*NumSamples,NumParams]);
    B = reshape(CholCovarCWLS' \ Z', NumSeries*NumSamples,1);
    Param = A\B;

    % Step 6 - combined E and CM step to estimate covariance parameters
    if celldesign && SingleDesign
        Mean = Design{1} * Param;
        Means = [];
    else
        Means = reshape(Xdesign*Param, [NumSeries, NumSamples]);
    end

    Covar0 = Covar;
    Covar = zeros(NumSeries,NumSeries);

    if SingleDesign
        Resid = Data - repmat(Mean',size(Z,1),1);
    else
        Resid = Data - Means';
    end
    CovAdj = zeros(NumSeries, NumSeries);

    if ~isempty(partialrows)  % always empty for mvn method
        WarnState = warning('off','MATLAB:nearlySingularMatrix');
        Z = zeros(1,NumSeries);
        for i = partialrows
            if ~SingleDesign
                Mean = Means(:,i);
            end

            [mX, mY, CXX, CXY, CYY] = ecmpart(Data(i,:), Mean, Covar0);

            if ~isempty(mY)
                CovAdj(:) = 0;
                if isempty(mX)
                    Z = Data(i,:);
                else
                    P = isnan(Data(i,:));
                    Q = ~P;
                    Y = Data(i,Q)';

                    Z(P) = mX + CXY * (CYY \ (Y - mY));
                    Z(Q) = Y;
                    CovAdj(P,P) = CXX - CXY * (CYY \ CXY');
                end

                Resid(i,:) = Z - Mean';

                Covar = Covar + CovAdj;
            end
        end
        warning(WarnState);
    end

    Covar = (Covar + Resid'*Resid) / Count;
    if isdiagonal
        Covar = diag(diag(Covar));
    end

    % Step 7 - evaluate objective and test for convergence
    %          update chol(covar) for non-lsq method only
    if EstMethod == 1  % cwls method
        Objective = statecmobj(Design,Data,Param,CovarCWLS,Resid,CholCovarCWLS,cwlsdiagonal);
    else
        [CholCovarCWLS, CholState] = chol(Covar);
        if CholState>0
            if EstMethod == 3
                error(message('stats:statmvnrobj:NonPosDefCov'));
            else
                error(message('stats:statecmobj:NonPosDefCov'));
            end
        end
        if EstMethod == 3     % mvn method
            Objective = statmvnrobj(Data,Design,Param,Covar,Resid,CholCovarCWLS,isdiagonal);
        else
            Objective = statecmobj(Design,Data,Param,Covar,Resid,CholCovarCWLS,isdiagonal);
        end
    end

    if Iter > 1
        TestObj = Objective - Objective0;
        TestParam = norm(Param - Param0)/sqrt(max(1,NumParams)); % force to zero if no params

        EpsObj = TolObj * (1 + abs(Objective));
        EpsParam = TolParam * (seps + norm(Param));

        if ((TestObj >= 0.0) && (TestObj < EpsObj)) && (TestParam < EpsParam)
            break
        end
    end

    Objective0 = Objective;
    Param0 = Param;
    
    if ~isempty(OutFun)
        str.Covar = Covar;
        str.iteration = Iter;
        str.fval = Objective;
        if OutFun(Param,str,'iter')
            disp(getString(message('stats:mvregress:TerminatedByUser')));
            return
        end
    end

    if (Iter == MaxIter) && (MaxIter > 1) && (TolObj > 0.0) && (TolParam > 0.0)
        warning(message('stats:mvregress:EarlyTermination'));
        break
    end
    
    if EstMethod == 3 && NumSeries == 1
        break
    end
end

if ~isempty(OutFun)
    str.Covar = Covar;
    str.iteration = Iter;
    str.fval = Objective;
    OutFun(Param,str,'done'); 
end

% Restore resids to proper size
if nargout>=3 && ~all(goodrows);
    Resid = resizeresids(Resid,goodrows);
end

% Compute parameter var/covar if requested
if nargout>=4
    % Fisher information matrix should be computed using the input covariance
    % for the 'cwls' method.  For other methods, set the cwls variables so that
    % the output covariance is used.
    if EstMethod ~= 1
         CovarCWLS = Covar;
         cwlsdiagonal = isdiagonal;
    end
    Info = statecmmvnrfish(Data,Design,CovarCWLS,VarType,VarFormat,CholCovarCWLS,cwlsdiagonal);
    VarParam = inv(Info);
    VarParam = (VarParam+VarParam')/2; % force to be exactly symmetric
end

if ~SingleDesign && ~iscell(Design0)
    Param = reshape(Param,[NumParams0 NumSeries]);
end

% -------------------------------------------------
function [mX, mY, CXX, CXY, CYY] = ecmpart(z, m, C)
%ECMPART Partitioning function for missing data algorithms
% Private routine to partition a mean vector m and covariance matrix C
% according to the pattern of NaNs (missing values) in an input vector z.

P = isnan(z);
Q = ~P;

mX = m(P);
mY = m(Q);

CXX = C(P,P);
CXY = C(P,Q);
CYY = C(Q,Q);

function outr = resizeresids(inr,goodrows)
outr = nan(length(goodrows),size(inr,2));
outr(goodrows,:) = inr;
