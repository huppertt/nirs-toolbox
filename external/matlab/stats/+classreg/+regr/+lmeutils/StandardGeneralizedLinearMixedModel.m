classdef StandardGeneralizedLinearMixedModel < classreg.regr.lmeutils.StandardLinearLikeMixedModel
%  This feature is intended for internal use only and is subject to change
%  at any time without warning.

%StandardGeneralizedLinearMixedModel - This is an internal utility to fit
%   Generalized Linear Mixed Effects (GLME) models in standard form. This
%   feature is intended for internal use only and is subject to change at
%   any time without warning.
%%
% $P(y_i | b) \sim Distr(\mu_i, \frac{\sigma^2}{w_i}, \phi)$
%
% $g(\mu) = X \beta + Z b + \delta$
% 
% * $y$ is a $N \times 1$ observed vector
% * $X$ is a $N \times p$ fixed effects design matrix
% * $\beta$ is a $p \times 1$ fixed effect vector
% * $Z$ is a $N \times q$ random effects design matrix
% * $b$ is a $q \times 1$ random effect vector
% * $\delta$ is a $N \times 1$ offset vector
% * $b \sim N(0, \sigma^2 D(\theta))$
% * $\Psi = \sigma^2 D(\theta)$
% * $Distr$ is a specified conditional distribution of $y | b$
% * $\mu$ is conditional mean of $y$ given $b$
% * $\sigma^2$ is the dispersion parameter
% * $\phi$ is a vector of additional parameters that appear in $Distr$
% * $g$ is the link function
% * $\phi$ is non-empty in some cases depending on $Distr$ (e.g., negative binomial distribution)
%
%   The purpose of this class is to:
%   (1) Fit GLME models in standard form by PL or ML.
%   (2) Refit the GLME model after changing some inputs.
%   (3) Compute the estimate betaHat of beta.
%   (4) Compute the estimate thetaHat of theta.
%   (5) Compute the estimate sigmaHat of sigma and phiHat of phi.
%   (6) Compute the empirical Bayes predictor (EBP) bHat of b.
%   (7) Estimate covariance of betaHat.
%   (8) Perform hypothesis tests on beta.
%   (9) CIs for contrasts of beta.
%  (10) Compute covariance of [beta;b] given y, thetaHat, sigmaHat, phiHat.
%  (11) Peform hypothesis tests on [beta;b].
%  (12) CIs for contrasts of [beta;b].
%  (13) Estimate covariance of [thetaHat;log(sigmaHat);phiHat].
%  (14) Estimate covariance of [etaHat;log(sigmaHat);phiHat] where etaHat 
%       is the Natural parameter vector for [thetaHat;log(sigmaHat)].
%  (15) Compute confidence intervals on the canonical parameter vector
%       [heta;sigma;phi].
%  (16) Compute conditional/marginal fitted values.
%  (17) Compute conditional/marginal residuals of various types.
%  (18) Generate random data from fitted model.
%  (19) Make predictions on new data.
%  (20) Compute the posterior covariance of b given y, betaHat, thetaHat, sigmaHat, phiHat.
%  (21) Compute various goodness of fit criterion.
%
%   StandardGeneralizedLinearMixedModel methods:
%      StandardGeneralizedLinearMixedModel - Create/fit a GLME model in standard form.
%      refit - Fit/refit a previously initialized GLME model.
%      initstats - Make a fitted GLME model ready for stats.
%      fitted - Get fitted values.
%      residuals - Convenience method to get various types of residuals.
%      modelCriterion - Compute various goodness of fit criterion.
%      fixedEffects - Summary stats on fixed effects.
%      randomEffects - Summary stats on random effects.
%      covarianceParameters - Summary stats on covariance parameters.
%      betaTTest - Perform T-test on beta.
%      betaFTest - Perform F-test on beta.
%      betaCI - Get confidence intervals on c'*beta.
%      betaBTTest - Perform T-test on [beta;b].
%      betaBFTest - Perform F-test on [beta;b].
%      betaBCI - Get confidence intervals on c'*[beta;b].
%      random - Generate random data from fitted model.
%      predict - Make predictions on new data.
%      predictTable - Summary table of predictions on new data.
%      postCovb - Posterior covariance of random effects.
%      hetaSigmaCI - Confidence intervals on covariance parameters (Canonical scale).
%      dfBetaTTest - Satterthwaite DF for T-test on beta.
%      dfBetaFTest - Satterthwaite DF for F-test on beta.
%      dfBetaBTTest - Satterthwaite DF for T-test on [beta;b].
%      dfBetaBFTest - Satterthwaite DF for F-test on [beta;b].
%      etaLogSigmaCI - Confidence intervals on natural parameters.
%      randomb - Generate random b vectors.
%      storeCovBetaHatBHat - Store covariance of [betaHat;bHat] in object.                                       
%      unstoreCovBetaHatBHat - Remove covariance of [betaHat;bHat] from object. 
%      turnOffOptimizerDisplay - Turn off optimizer display.
%
%   StandardGeneralizedLinearMixedModel properties:
%      y - Response vector.
%      X - Fixed effects design matrix
%      Z - Random effects design matrix.
%      Psi - An object of type CovarianceMatrix encapsulating $\sigma^2 D$.
%      FitMethod - Method used to fit the GLME.  
%      Distribution - Name for the conditional distribution of y given b.
%      Link - A structure containing link function info.
%      Offset - Offset vector.
%      BinomialSize - Trials per observation for binomial response.
%      PriorWeights - Vector of prior observation weights.
%      DispersionFixed - True if dispersion parameter is fixed at 1.0.
%      VarianceFunction - A structure containing variance function info.
%      N - Number of rows in y, X, Z and delta.
%      p - Number of columns in X.       
%      q - Number of columns in Z.
%      rankX - Rank of the fixed effects design matrix X.
%      Optimizer - Selected optimizer.
%      OptimizerOptions - Options for the selected optimizer.
%      CheckHessian - Flag for Hessian checks at convergence.
%      PLIterations - Maximum number of PL iterations.
%      PLTolerance - Relative tolerance for PL iterations.
%      MuStart - Conditional mean of y given b used to initialize PL iterations.
%      MuBound - Bounds imposed on the conditional mean of y given b during iterations.
%      InitPLIterations - Initial number of PL iterations before starting ML estimation.
%      EBMethod - Name of method to use for EB sub optimization.
%      EBOptions - Options for EB sub optimization.
%      CovarianceMethod - Method for computing covariance of parameter estimates.
%      UseSequentialFitting - True for sequential fitting strategy for ML estimation.
%      InitializationMethod - Method used to initialize the first PL iteration.
%      slme - A StandardLinearMixedModel object from final PL iteration.
%      betaHat - Estimated fixed effects vector.        
%      bHat - Estimated EBP of random effects vector b.
%      DeltabHat - Estimated normalized EBP of random effects vector b.
%      sigmaHat - Estimated square root of dispersion parameter.       
%      thetaHat - Estimated unconstrained parameter for matrix D.
%      phiHat - Estimated additional parameter in the distribution of P(y|b).
%      loglikHat - Maximized log likelihood or restricted log likelihood.   
%      loglikHatPseudoData - Maximized log likelihood or restricted log likelihood of pseudo data from final PL iteration.
%      covbetaHat - Estimated covariance of betaHat.
%      covthetaHatlogsigmaHat - Estimated covariance on unconstrained scale of [thetaHat;log(sigmaHat)] or [thetaHat;log(sigmaHat);phiHat].       
%      covetaHatlogsigmaHat - Estimated covariance on Natural scale of [etaHat;log(sigmaHat)] or [etaHat;log(sigmaHat);phiHat].        
%      covbetaHatbHat - Estimated covariance of [beta;b] given y, thetaHat, sigmaHat, phiHat.
%      isSigmaFixed - True if residual standard deviation sigma is fixed.
%      sigmaFixed - The specified fixed value of sigma.
%      isFitToData - Flag to mark a fitted StandardGeneralizedLinearMixedModel object.
%      isReadyForStats - Flag to mark the object ready for doing stats.
%      UseAMDPreordering - Logical flag indicating whether to use AMD preordering during model fitting.
%      AMDOrder - A q-by-1 integer vector containing the precomputed AMD order.
%      NewtonStepMethod - Method used to compute the Newton step during posterior mode estimation.    

%   Copyright 2013-2014 The MathWorks, Inc.    

% Public properties supplied by the user.
    properties (GetAccess=public, SetAccess=public)
%y - N-by-1 response vector used to fit the GLME.        
        y

%X - N-by-p fixed effects design matrix used to fit the GLME.        
        X

%Z - N-by-q random effects design matrix used to fit the GLME.       
        Z
        
%Psi - An object of type CovarianceMatrix representing the matrix Psi =
%   sigma^2 * D(theta). When optimization completes, Psi will contain the
%   optimal theta and sigma.
        Psi
        
%FitMethod - The method used to fit the GLME model:
%                  'MPL' - Maximum pseudo likelihood
%                'REMPL' - Restricted maximum pseudo likelihood
%   'ApproximateLaplace' - Maximum likelihood using approximate Laplace
%                          approximation with fixed effects profiled out
%              'Laplace' - Maximum likelihood using Laplace approximation
%           'Quadrature' - Maximum likelihood using adaptive Gauss-Hermite 
%                          quadrature (TODO)
        FitMethod                
    end
    
% Other read only public properties.    
    properties (GetAccess=public, SetAccess=protected)
%Distribution - A string containing the name of the conditional
%   distribution of y given b.
        Distribution        
        
%Link - A structure containing g, its inverse, 1st and 2nd derivatives:
%        Name                 Name of the link function.
%        Link                 The function that defines g.
%        Derivative           Derivative of g.
%        SecondDerivative     Second derivative of g.
%        Inverse              Inverse of g.
        Link
        
%Offset - N-by-1 offset vector used to fit the GLME.
        Offset
        
%BinomialSize - N-by-1 vector of integers containing the number of trials
%   using which the proportions in y were computed.
        BinomialSize
        
%PriorWeights - N-by-1 vector of prior observation weights. For binomial 
%   and Poisson distributions, these are positive integers.
        PriorWeights        
        
%DispersionFixed - true if the dispersion parameter needs to be fixed at 
%   1.0 and false otherwise. This option only makes sense for PL based
%   methods and only applies to binomial and Poisson distributions.
        DispersionFixed
        
%VarianceFunction - Structure with variance function and its derivative:
%        VarianceFunction     Variance function for y given b as a function
%                             of conditional mean mu of y given b.
%        Derivative           Derivative of variance function for y given b
%                             as a function of conditional mean mu of y 
%                             given b.
        VarianceFunction
    end
    
% Public read only properties storing size information.   
    properties (GetAccess=public, SetAccess=protected)
%N - Number of observations. This is the number of rows in y, X and Z.
        N
        
%p - Number of fixed effects. This is the number of columns in X. 
        p
        
%q - Number of random effects. This is the number of columns in Z.
        q        
        
%rankX - Rank of the fixed effects design matrix.
        rankX
    end
  
% Public read only properties storing optimization related info. 
    properties(GetAccess=public, SetAccess=protected)
%Optimizer - A string containing the name of the optimizer.        
        Optimizer = 'fminsearch';
        
%OptimizerOptions - A structure containing the optimizer options.        
        OptimizerOptions = struct([]);        

%CheckHessian - A logical scalar indicating whether to perform the Hessian 
%   checks in initstats. Hessian checks can be expensive for GLME models 
%   and so this functionality is currently not available. (TODO)
        CheckHessian = false;        
        
%PLIterations - A positive integer specifying the maximum number of pseudo 
%   likelihood (PL) iterations. Applies only if 'FitMethod' is 'MPL' or 
%   'REMPL'.
        PLIterations     
        
%PLTolerance - A numeric, real scalar indicating the relative convergence
%   tolerance for terminating PL iterations.
        PLTolerance        
        
%MuStart - A N-by-1 vector providing a starting value for the conditional 
%   mean of y given b to initialize PL iterations. Legal values of elements 
%   in MuStart are as follows:
%
%                 Distribution          Legal values
%                 'normal'              (-Inf,Inf)
%                 'binomial'            (0,1)
%                 'poisson'             (0,Inf)
%                 'gamma'               (0,Inf)
%                 'inverse gaussian'    (0,Inf)        
%
%   If MuStart is empty, a default initialization for MuStart is used.
        MuStart

%MuBound - A structure with two fields TINY and BIG that control the 
%   distribution specific checks each time the conditional mean mu is
%   calculated from the linear predictor eta by appling the inverse link.
%   Elements of mu are required to take on legal values as defined below:
%
%                 Distribution          Legal values
%                 'normal'              (-BIG,BIG)
%                 'binomial'            (TINY,1-TINY)
%                 'poisson'             (TINY,BIG)
%                 'gamma'               (TINY,BIG)
%                 'inverse gaussian'    (TINY,BIG)  
%
%   where TINY = MuBound.TINY and BIG = MuBound.BIG. Default value of TINY
%   is eps and BIG is Inf. When mu does not satisy the above constraints,
%   mu is redefined by setting its elements smaller than max(eps,TINY) to
%   max(eps,TINY) and elements bigger than min(realmax,BIG) to
%   min(realmax,BIG).
        MuBound = struct('TINY',eps,'BIG',Inf);        
        
%InitPLIterations - Initial number of PL iterations used to initialize ML 
%   based methods like 'Laplace', 'ApproximateLaplace' and 'Quadrature'.
        InitPLIterations        
        
%EBMethod - A string containing the optimization method to use for
%   estimating the empirical Bayes (EB) predictors of random effects.
%   Choices are 'Default', 'LineSearchNewton', 'LineSearchModifiedNewton'
%   and 'TrustRegion2D'. The 'Default' 'EBMethod' is very similar to
%   'LineSearchNewton' but uses a different convergence criterion and does
%   not display iterative progress. 'Default' and 'LineSearchNewton' may
%   fail for non-canonical link functions. In these cases, 'TrustRegion2D'
%   is the recommended 'EBMethod'.
        EBMethod     
        
%EBOptions - A structure created using statset to control the EB
%   optimization. EBOptions has the following fields:
%
%                 'TolFun'      Relative tolerance on the gradient norm.
%                               Default is 1e-6.
%                 'TolX'        Absolute tolerance on the step size.
%                               Default is 1e-8.
%                 'MaxIter'     Maximum number of iterations. Default is
%                               100.
%                 'Display'     'off', 'iter' or 'final'. Default is 'off'.
%
%   For 'EBMethod' equal to 'Default', 'TolFun' is the relative tolerance
%   on the linear predictor of the model and the option 'Display' does not
%   apply.
        EBOptions                          
               
%CovarianceMethod - A string indicating the method used to compute the 
%   covariance matrix of estimated parameters. Choices are 'Conditional'
%   and 'JointHessian'. The 'Conditional' method is fast and computes an
%   approximate covariance of fixed effects conditional on the estimated
%   covariance parameters equal to the true values. The 'JointHessian'
%   method computes the joint Hessian of the Laplacian log likelihood with
%   respect to all model parameters and uses this to compute the covariance
%   of fixed effects and covariance parameters. Default is 'Conditional'.
%   For the 'Conditional' method, the covariance of covariance parameters
%   is not computed.
        CovarianceMethod
        
        
%UseSequentialFitting - Either true or false. If false, all ML methods are 
%                       initialized using 1 or more PL iterations. If true,
%                       the initial values from PL iterations are refined
%                       with 'ApproximateLaplace' for 'Laplace' based
%                       fitting and with 'ApproximateLaplace' followed by
%                       'Laplace' for 'Quadrature' based fitting.
        UseSequentialFitting
    end
    
% Public property storing initialization info.    
    properties(GetAccess=public, SetAccess=public)
%InitializationMethod - A string containing the method used to compute the 
%   initial value of theta in the first PL iteration.        
        InitializationMethod = 'default';        
    end
    
% Public read only property applicable to PL based methods.    
    properties(GetAccess=public, SetAccess=protected)        
%slme - A StandardLinearMixedModel object containing the final LME fit if 
%   'FitMethod' is 'MPL' or 'REMPL' and empty otherwise.        
        slme        
    end
    
% Public read only properties for estimated GLME model parameters.   
    properties (GetAccess=public, SetAccess=protected)
%betaHat - Estimated fixed effects vector.
        betaHat
        
%bHat - Estimated EBP of random effects vector.
        bHat
        
%DeltabHat - Estimated normalized EBP of random effects vector.
        DeltabHat        
        
%sigmaHat - Estimated square root of dispersion parameter.
        sigmaHat
        
%thetaHat - Estimated vector of unconstrained parameters for matrix D.
        thetaHat        
        
%phiHat - Estimated vector of extra parameters in the conditional distribution of y given b (beyond dispersion).        
        phiHat

%loglikHat - Maximized log likelihood or maximized restricted log likelihood.       
        loglikHat    
                
%loglikHatPseudoData - Maximized log likelihood or maximized restricted log 
%   likelihood of pseudo data using the final LME model if 'FitMethod' is
%   'MPL' or 'REMPL'.
        loglikHatPseudoData      
    end
        
% Public read only properties for stats info on GLME model parameters.
    properties (GetAccess=public, SetAccess=protected)        
%covbetaHat - Estimated covariance of betaHat.        
        covbetaHat
                
%covthetaHatlogsigmaHat - Estimated covariance of [thetaHat;log(sigmaHat)] or [thetaHat;log(sigmaHat);phiHat]. 
        covthetaHatlogsigmaHat      
        
%covetaHatlogsigmaHat - Estimated covariance of [etaHat;log(sigmaHat)] or [etaHat;log(sigmaHat);phiHat].
%   etaHat is the Natural parameter vector corresponding to thetaHat.        
        covetaHatlogsigmaHat
        
%covbetaHatbHat - Estimated covariance of [beta;b] given y, thetaHat, sigmaHat, phiHat.
        covbetaHatbHat
    end
        
% Public read only properties storing whether sigma is fixed and if so its value.
    properties (GetAccess=public, SetAccess=protected)
%isSigmaFixed - True if square root of dispersion parameter is fixed.
        isSigmaFixed = false;
        
%sigmaFixed - Scalar fixed value for square root of dispersion parameter.
        sigmaFixed = NaN;
    end    
    
% Public read only properties storing state information.
    properties (GetAccess=public, SetAccess=protected)
%isFitToData - True if this is a fitted StandardGeneralizedLinearMixedModel object.
        isFitToData = false;
        
%isReadyForStats - True if this StandardGeneralizedLinearMixedModel object is ready for stats.
        isReadyForStats = false;
    end
    
% Public read only properties related to AMD ordering.    
    properties (GetAccess=public, SetAccess=protected)                
%UseAMDPreordering - A logical flag indicating whether AMD preordering 
%   should be used during model fitting. This is applicable only to ML
%   based methods. If true, an AMD preordering of U'*C*U + I_q is computed
%   and stored in the object. This preordering is then applied before
%   taking the Cholesky factorization of U'*C*U + I_q.
        UseAMDPreordering
        
%AMDOrder - A q-by-1 integer vector storing the AMD preordering for use in 
%   posterior mode estimation.        
        AMDOrder
        
%NewtonStepMethod - A string indicating the method to use for computing 
%   the Newton step during posterior mode estimation. Valid values are:
%
%       'Cholesky'  - Attempt Cholesky factorization of (U'*C*U + I_q) 
%                     first (possibly with AMD preordering) and if not 
%                     successful, use backslash \.
%       'Backslash' - Always use backslash \.
        NewtonStepMethod
    end

% Private properties related to AMD ordering.    
    properties (Access=private)
%NewtonStepMethodCode - An integer indicating the method to use for 
%   computing the Newton step during posterior mode estimation. The mapping
%   between 'NewtonStepMethod' and 'NewtonStepMethodCode' is as follows:
%
%       'NewtonStepMethod'          'NewtonStepMethodCode'
%          'Cholesky'                        1
%          'Backslash'                       2
        NewtonStepMethodCode        
    end

% Private property related to detection of badly scaled weights in PL.    
    properties (Access=private)
%HaveWarnedAboutBadlyScaledPLWeights - Logical flag indicating whether we
%   have already warned about badly scaled weights during PL iterations or
%   not. This is useful for displaying a warning related to badly scaled PL
%   weights once instead of at every PL iteration.
        HaveWarnedAboutBadlyScaledPLWeights = false;
    end
    
% Private property related to PL optimizer display.
    properties (Access=private)
%ShowPLOptimizerDisplay - Logical flag indicating whether we want to show 
%   the iterative display from PL optimizer.    
        ShowPLOptimizerDisplay = false;        
    end
    
% Internal constants.
    properties (Constant=true, Hidden=true)       
%AllowedFitMethods - The fitting methods that we currently support.        
        AllowedFitMethods = {'mpl','rempl','approximatelaplace','laplace','quadrature'};        
        
%AllowedDistributions - A list of allowed conditional distributions for y given b.
        AllowedDistributions = {'normal','gaussian','binomial','poisson','gamma','inversegaussian','inverse gaussian'};
        
%AllowedStandardLinks - A list of allowed 'standard' link functions.
        AllowedLinks = {'identity','log','logit','probit','comploglog','loglog','reciprocal'};     
        
%AllowedEBMethods - A list of allowed optimization methods for estimating
%   empirical Bayes predictors of random effects.
        AllowedEBMethods = {'Default','LineSearchNewton','LineSearchModifiedNewton','TrustRegion2D','fsolve'}; 
               
%AllowedCovarianceMethods - A list of allowed methods for computing
%   covariance of estimated parameters.
        AllowedCovarianceMethods = {'Conditional','JointHessian'};
        
%AllowedNewtonStepMethods - A list of allowed methods for computing the 
%   Newton step during posterior mode estimation.
        AllowedNewtonStepMethods = {'Cholesky','Backslash'};
    end                
    
% Set methods for X, y, Z, Psi and FitMethod.   
    methods
        
        function sglme = set.y(sglme,newy)
        
            if ~isempty(sglme.y)
                % (1) Validate newy.
                newy = validatey(sglme,newy);
                
                % (2) Invalidate the fit.
                sglme = invalidateFit(sglme);
            end
            
            % (3) Set the new y.
            sglme.y = newy;
            
        end % end of set.y
        
        function sglme = set.X(sglme,newX)
            
            if ~isempty(sglme.X)
                % (1) Validate newX.
                newX = validateX(sglme,newX);
                
                % (2) Invalidate the fit.
                sglme = invalidateFit(sglme);
            end
            
            % (3) Set the new X.
            sglme.X = newX;
            
        end % end of set.X
        
        function sglme = set.Z(sglme,newZ)
            
            if ~isempty(sglme.Z)
                % (1) Validate newZ.
                newZ = validateZ(sglme,newZ);
                
                % (2) Invalidate the fit.
                sglme = invalidateFit(sglme);
            end            
            
            % (3) Set the new Z.
            sglme.Z = newZ;
            
        end % end of set.Z
        
        function sglme = set.Psi(sglme,newPsi)
            
            if ~isempty(sglme.Psi)
                % (1) Validate newPsi.
                newPsi = validatePsi(sglme,newPsi);
                
                % (2) Invalidate the fit.
                sglme = invalidateFit(sglme);
            end
            
            % (3) Set the new Psi.
            sglme.Psi = newPsi;
            
        end % end of set.Psi
        
        function sglme = set.FitMethod(sglme,newFitMethod)
            
            if ~isempty(sglme.FitMethod)
                % (1) Validate newFitMethod.
                newFitMethod = validateFitMethod(sglme,newFitMethod);
                
                % (2) Invalidate the fit.
                sglme = invalidateFit(sglme);
            end
            
            % (3) Set the new FitMethod.
            sglme.FitMethod = newFitMethod;            
            
        end % end of set.FitMethod
        
    end
    
% Protected validation/utility methods.
    methods (Access=protected)  
        
        function FitMethod = validateFitMethod(sglme,FitMethod)
            
            % (1) FitMethod is a string from the list AllowedFitMethods.
            FitMethod = internal.stats.getParamVal(FitMethod,sglme.AllowedFitMethods,'FitMethod');
            
        end % end of validateFitMethod.                                       
                               
        function offset = validateOffset(sglme,offset,N) %#ok<INUSL>
%offset = validateOffset(sglme,offset,N) takes a vector offset and an 
%   integer N, and returns the validated offset vector. offset must be a
%   numeric, real vector of size N-by-1.

            % (1) Ensure that N is a scalar integer.
            assert(internal.stats.isScalarInt(N));
            
            % (2) Ensure that offset is a N-by-1 numeric, real vector.
            isok = isnumeric(offset) & isreal(offset) & ...
                isvector(offset) & size(offset,1) == N;            
            if ~isok
                error(message('stats:classreg:regr:lmeutils:StandardGeneralizedLinearMixedModel:BadOffset',N));
            end            
            
        end % end of validateOffset.        
        
        function binomialsize = validateBinomialSize(sglme,binomialsize,N) %#ok<INUSL>
%binomialsize = validateBinomialSize(sglme,binomialsize,N) takes a vector 
%   binomialsize and an integer N and returns the validated binomialsize.
%   binomialsize must be a vector of integers > 0 and it must be of size
%   N-by-1.
           
            % (1) Ensure that N is a scalar integer.
            assert(internal.stats.isScalarInt(N));
            
            % (2) Ensure that binomialsize is a N-by-1 vector of positive
            % integers.
            isok = internal.stats.isIntegerVals(binomialsize,1) & ...
                isvector(binomialsize) & size(binomialsize,1) == N;
            if ~isok
                error(message('stats:classreg:regr:lmeutils:StandardGeneralizedLinearMixedModel:BadBinomialSize',N));
            end                                
            
        end % end of validateBinomialSize.                                      
        
        function weights = validateWeights(sglme,weights,N,distribution) %#ok<INUSL>
%weights = validateWeights(sglme,weights,N,distribution) takes a vector of 
%   observation weights, an integer N and a string distribution. weights is
%   validated and returned. We expect weights to be of size N-by-1 with
%   positive integer elements if distribution is 'binomial' or 'poisson' 
%   and a numeric real vector of size N-by-1 otherwise.

            % (1) Ensure that N is a scalar integer.
            assert(internal.stats.isScalarInt(N));

            % (2) Ensure that distribution is a string.
            assert(internal.stats.isString(distribution));
            
            % (3) Ensure that weights is a numeric, real vector of size
            % N-by-1. If distribution is 'binomial' or 'poisson' then
            % weights must have positive integer values.            
            if any(strcmpi(distribution,{'binomial','poisson'}))
                isok = internal.stats.isIntegerVals(weights,1) & ...
                    isvector(weights) & size(weights,1) == N;
            else
                isok = isnumeric(weights) & isreal(weights) ...
                    & isvector(weights) & size(weights,1) == N;
            end
            if ~isok
                error(message('stats:classreg:regr:lmeutils:StandardGeneralizedLinearMixedModel:BadWeights',N,N));
            end                        
            
        end % end of validateWeights.
                
        function validateyRange(sglme,y,binomialsize,weights,distribution)
%validateyRange(sglme,y,binomialsize,weights,distribution) takes a N by 1 
%   response vector y, a N by 1 vector binomialsize, a N by 1 vector
%   weights, a string distribution and validates that values in y are
%   sensible given the selected distribution, binomialsize (if applicable)
%   and weights. Here's what is checked:
%
%   Distribution                      Constraint on y
%   Binomial                          y in [0,1], weights.*binomialsize.*y is an integer
%   Poisson                           y >= 0, weights.*y is an integer
%   Gamma                             y > 0
%   Inverse Gaussian                  y > 0
%   Normal                            y in (-Inf,Inf)

            % (1) Ensure that y is of size N-by-1, numeric, real vector.
            assert(iscolumn(y) & isnumeric(y) & isreal(y));
            N = size(y,1);
            
            % (2) Ensure that distribution is a valid distribution.
            distribution = internal.stats.getParamVal(distribution,sglme.AllowedDistributions,'Distribution'); 
            
            % (3) Ensure that binomialsize is valid.
            binomialsize = validateBinomialSize(sglme,binomialsize,N);
            
            % (4) Ensure that weights are valid.
            weights = validateWeights(sglme,weights,N,distribution);
            
            % (5) Check the values in y.
            switch lower(distribution)              
                case 'binomial'                    
                    counts = weights.*binomialsize.*y;                                                            
                    isok = all(y >= 0 & y <= 1) & max(abs(counts - round(counts))) <= sqrt(eps);
                    if ~isok
                        error(message('stats:classreg:regr:lmeutils:StandardGeneralizedLinearMixedModel:BadBinomialY'));
                    end                    
                case 'poisson'
                    counts = weights.*y;
                    isok = all(y >= 0) & max(abs(counts - round(counts))) <= sqrt(eps);
                    if ~isok
                        error(message('stats:classreg:regr:lmeutils:StandardGeneralizedLinearMixedModel:BadPoissonY'));
                    end
                case 'gamma'
                    isok = all(y > 0);
                    if ~isok
                        error(message('stats:classreg:regr:lmeutils:StandardGeneralizedLinearMixedModel:BadGammaY'));
                    end                    
                case {'inverse gaussian','inversegaussian'}
                    isok = all(y > 0);
                    if ~isok
                        error(message('stats:classreg:regr:lmeutils:StandardGeneralizedLinearMixedModel:BadInverseGaussianY'));
                    end                    
                case {'normal','gaussian'}                                        
            end            
            
        end % end of validateyRange.
        
        function linkSpec = defaultLink(sglme,distribution)
%linkSpec = defaultLink(sglme,distribution) takes a string distribution 
%   that represents one of the allowed distributions and returns the
%   default (canonical) link for the distribution.

            % (1) Ensure that distribution is legal.
            distribution = internal.stats.getParamVal(distribution,sglme.AllowedDistributions,'Distribution');      

            % (2) Get the default link.
            switch lower(distribution)            
                case {'normal','gaussian'}
                    linkSpec = 'identity';
                case 'binomial'
                    linkSpec = 'logit';
                case 'poisson'
                    linkSpec = 'log';
                case 'gamma'
                    linkSpec = -1;
                case {'inverse gaussian','inversegaussian'}
                    linkSpec = -2;
            end
                
        end % end of defaultLink.
        
        function varStruct = varianceFunction(sglme,distribution)
%varStruct = varianceFunction(sglme,distribution) takes a string 
%   distribution that is one of the supported distributions and returns a
%   struct varStruct with the following two fields:
%
%        VarianceFunction     Variance function for y given b as a function
%                             of conditional mean mu of y given b.
%        Derivative           Derivative of variance function for y given b
%                             as a function of conditional mean mu of y 
%                             given b.
%
%   The purpose of this helper method is to initialize the field
%   VarianceFunction.

            % (1) Ensure that distribution is legal.
            distribution = internal.stats.getParamVal(distribution,sglme.AllowedDistributions,'Distribution');      

            % (2) Get the variance function info.
            switch lower(distribution)            
                case {'normal','gaussian'}
                    varStruct.VarianceFunction = @(mu) ones(size(mu));
                    varStruct.Derivative = @(mu) zeros(size(mu));
                case 'binomial'
                    varStruct.VarianceFunction = @(mu) mu.*(1-mu);
                    varStruct.Derivative = @(mu) 1-2*mu;
                case 'poisson'
                    varStruct.VarianceFunction = @(mu) mu;
                    varStruct.Derivative = @(mu) ones(size(mu));
                case 'gamma'
                    varStruct.VarianceFunction = @(mu) mu.^2;
                    varStruct.Derivative = @(mu) 2*mu;
                case {'inverse gaussian','inversegaussian'}
                    varStruct.VarianceFunction = @(mu) mu.^3;
                    varStruct.Derivative = @(mu) 3*(mu.^2);
            end
            
        end % end of varianceFunction.
                
        function pliterations = validatePLIterations(sglme,pliterations) %#ok<INUSL>
%pliterations = validatePLIterations(sglme,pliterations) takes a potential
%   value for the number of PL iterations and returns the validated value.

            % (1) pliterations must be a positive scalar integer.
            isok = isscalar(pliterations) & ...
                internal.stats.isIntegerVals(pliterations,1);
            if ~isok
                error(message('stats:classreg:regr:lmeutils:StandardGeneralizedLinearMixedModel:BadPLIterations'));
            end
            
        end % end of validatePLIterations.
        
        function pltolerance = validatePLTolerance(sglme,pltolerance) %#ok<INUSL>
%pltolerance = validatePLTolerance(sglme,pltolerance) takes a potential
%   value for the tolerance for PL iterations and returns the validated
%   value.
        
            % (1) pltolerance must be a numeric, real scalar.
            isok = isscalar(pltolerance) & ...
                isnumeric(pltolerance) & isreal(pltolerance);
            if ~isok
                error(message('stats:classreg:regr:lmeutils:StandardGeneralizedLinearMixedModel:BadPLTolerance'));
            end
            
        end % end of validatePLTolerance.
        
        function mustart = validateMuStart(sglme,mustart,distribution,N)
%mustart = validateMuStart(sglme,mustart,distribution,N) takes a potential 
%   value of sglme.MuStart and validates it. Legal values of elements in 
%   mustart depend on the value of distribution as follows:
%
%                 Distribution          Legal values
%                 'normal'              (-Inf,Inf)
%                 'binomial'            (0,1)
%                 'poisson'             (0,Inf)
%                 'gamma'               (0,Inf)
%                 'inverse gaussian'    (0,Inf)      
%
%   mustart is also required to be of size N-by-1. Empty values of mustart
%   are acceptable.

            % (1) Ensure that N is sensible.
            assert(internal.stats.isIntegerVals(N,1) & isscalar(N));

            % (2) Ensure that distribution is sensible.
            assert(any(strcmpi(distribution,sglme.AllowedDistributions)));
            
            % (3) Check mustart.
            if ~isempty(mustart)
                sizeok = isvector(mustart) & (size(mustart,1) == N);
                switch lower(distribution)
                    case {'binomial'}
                        isok = sizeok & all(mustart > 0 & mustart < 1);                        
                    case {'poisson'}
                        isok = sizeok & all(mustart > 0 & mustart < Inf);
                    case {'gamma'}
                        isok = sizeok & all(mustart > 0 & mustart < Inf);
                    case {'inverse gaussian','inversegaussian'}
                        isok = sizeok & all(mustart > 0 & mustart < Inf);
                    case {'normal','gaussian'}
                        isok = sizeok & all(mustart > -Inf & mustart < Inf);                        
                end
                if ~isok
                    error(message('stats:classreg:regr:lmeutils:StandardGeneralizedLinearMixedModel:BadMuStart',N));
                end
            end                    
            
        end % end of validateMuStart.
        
        function dispersionfixed = setDispersionFixed(sglme,dispersionflag,distribution,fitmethod)
%dispersionfixed 
%        = setDispersionFixed(sglme,dispersionflag,distribution,fitmethod)
%   takes a logical scalar dispersionflag, a string distribution and a
%   string fitmethod and returns a logical scalar dispersionfixed such
%   that:
%
%   dispersionfixed                   Meaning
%       true                    GLME dispersion parameter is fixed at 1.0
%       false                   GLME dispersion parameter is NOT fixed
%
%   Here are the possible values of dispersionfixed depending on fitmethod
%   and distribution:
%
%                                fitmethod = {'MPL','REMPL'}      {'Laplace','ApproximateLaplace','Quadrature'}
%
%   distribution = 'binomial'          true/false                    true
%   distribution = 'poisson'           true/false                    true
%   Other values of distribution        false                        false
%
%   For fitmethod = {'MPL','REMPL'} we can force dispersionfixed = false if
%   dispersionflag = true.

            % (1) Basic checks on input.
            assert(isscalar(dispersionflag) & islogical(dispersionflag));
            assert(any(strcmpi(distribution,sglme.AllowedDistributions)));
            assert(any(strcmpi(fitmethod,sglme.AllowedFitMethods)));

            % (2) Fill default value into dispersionfixed.
            switch lower(distribution)
                case {'binomial','poisson'}
                    dispersionfixed = true;
                otherwise
                    dispersionfixed = false;
            end
            
            % (3) If dispersionflag is true and we have a binomial/poisson
            % distribution and if fitmethod is 'mpl' or 'rempl' then set
            % dispersionfixed to false.
            estdisp = (dispersionflag == true) & ...
                any(strcmpi(distribution,{'binomial','poisson'})) ...
                & any(strcmpi(fitmethod,{'mpl','rempl'}));
            if estdisp == true
                dispersionfixed = false;
            end            
            
        end % end of setDispersionFixed.        
                                
        function [ebmethod,eboptions] = validateEBParameters(sglme,ebmethod,eboptions,dfltEBOptions)
%[ebmethod,eboptions] = validateEBParameters(sglme,ebmethod,eboptions,dfltEBOptions)
%   takes a string ebmethod, a structure eboptions and a default structure
%   dfltEBOptions and validates ebmethod and eboptions. ebmethod must be a
%   string defined in sglme.AllowedEBMethods and eboptions must be a
%   structure.

            % (1) Validate ebmethod.
            ebmethod = internal.stats.getParamVal(ebmethod,sglme.AllowedEBMethods,'EBMethod');
                                    
            % (2) eboptions is a struct or 'optim.options.Fsolve' object.
            if ~isstruct(eboptions) && ~isa(eboptions,'optim.options.Fsolve')               
                error(message('stats:classreg:regr:lmeutils:StandardGeneralizedLinearMixedModel:BadEBOptions'));
            end
            
            % (3) There are 4 possible combinations of ebmethod/eboptions.
            if isstruct(eboptions) && strcmpi(ebmethod,'fsolve')
                % Filled out 'optim.options.Fsolve' object in eboptions.
                eboptions = sglme.convertOptionsToFSolveOptions(eboptions,dfltEBOptions);                                 
            elseif isstruct(eboptions) && ~strcmpi(ebmethod,'fsolve')
                % Filled out options struct in eboptions.
                eboptions = statset(dfltEBOptions,eboptions);
            elseif isa(eboptions,'optim.options.Fsolve') && strcmpi(ebmethod,'fsolve')
                % OK.
            elseif isa(eboptions,'optim.options.Fsolve') && ~strcmpi(ebmethod,'fsolve')
                % Can't supply 'optim.options.Fsolve' object for 
                % non-fsolve methods.
                error(message('stats:classreg:regr:lmeutils:StandardGeneralizedLinearMixedModel:BadEBOptions'));
            else
                % At least one of the 4 options must match.
                error(message('stats:classreg:regr:lmeutils:StandardGeneralizedLinearMixedModel:BadEBOptions'));
            end                                                 
            
        end % end of validateEBParameters.
        
        function initpliterations = validateInitPLIterations(sglme,initpliterations) %#ok<INUSL>
%initpliterations = validateInitPLIterations(sglme,initpliterations) takes
%   a potential value of initial number of PL iterations and validates it.

            % (1) initpliterations must be a scalar positive integer.
            isok = isscalar(initpliterations) & ...
                internal.stats.isIntegerVals(initpliterations,1);
            if ~isok
                error(message('stats:classreg:regr:lmeutils:StandardGeneralizedLinearMixedModel:BadInitPLIterations'));
            end

        end % end of validateInitPLIterations.
        
        function [mulowerbound,muupperbound] = validateMuBounds(sglme,mulowerbound,muupperbound) %#ok<INUSL>
%[mulowerbound,muupperbound] = validateMuBounds(sglme,mulowerbound,muupperbound)
%   takes potential values of mulowerbound and muupperbound and validates
%   them. mulowerbound must be a real, scalar in [0,1) and muupperbound
%   must be a real, positive scalar.

            % 1. Validate mulowerbound.
            ok = isnumeric(mulowerbound) & isreal(mulowerbound) & isscalar(mulowerbound);
            ok = ok & (mulowerbound >= 0 & mulowerbound < 1);
            if ~ok
                error(message('stats:classreg:regr:lmeutils:StandardGeneralizedLinearMixedModel:BadMuLowerBound'));
            end
            
            % 2. Validate muupperbound.
            ok = isnumeric(muupperbound) & isreal(muupperbound) & isscalar(muupperbound);
            ok = ok & (muupperbound > 0);
            if ~ok
                error(message('stats:classreg:regr:lmeutils:StandardGeneralizedLinearMixedModel:BadMuUpperBound'));
            end
            
            % 3. mulowerbound must be smaller than muupperbound.
            ok = mulowerbound < muupperbound;
            if ~ok
                error(message('stats:classreg:regr:lmeutils:StandardGeneralizedLinearMixedModel:BadMuLowerBound'));
            end

        end % end of validateMuBounds.
        
        function checkDistributionLinkCombination(sglme,distribution,linkStruct) %#ok<INUSL>
% Takes a validated string distribution that is one of 'normal','gaussian',
% 'binomial','poisson','gamma','inversegaussian','inverse gaussian' and a 
% validated linkStruct with the fields:
%
%        Name                 Name of the link function.
%        Link                 The function that defines g.
%        Derivative           Derivative of g.
%        SecondDerivative     Second derivative of g.
%        Inverse              Inverse of g.            
%
% and checks if the specified distribution and link combination is
% sensible. If not, a warning message is issued. The linkStruct.Name field
% can have these values: 'identity','log','logit','probit','comploglog',
% 'loglog','reciprocal','power',''

            % 1. The GLME linear predictor is unconstrained but the 
            % 'reciprocal' and 'power' links require the linear predictor
            % to be non-negative. 
            linkname = linkStruct.Name;
            if any(strcmpi(linkname,{'reciprocal','power'}))
                warning(message('stats:classreg:regr:lmeutils:StandardGeneralizedLinearMixedModel:BadDistLinkCombination1'));
            end

            % 2. Can't use the identity and log links for a binomial
            % distribution.
            if strcmpi(distribution,'binomial') && any(strcmpi(linkname,{'identity','log'}))
                warning(message('stats:classreg:regr:lmeutils:StandardGeneralizedLinearMixedModel:BadDistLinkCombination2'));
            end
            
            % 3. Can't use the identity link for poisson, gamma and inverse
            % gaussian distributions.
            if any(strcmpi(distribution,{'poisson','gamma','inversegaussian','inverse gaussian'})) && strcmpi(linkname,'identity')
                warning(message('stats:classreg:regr:lmeutils:StandardGeneralizedLinearMixedModel:BadDistLinkCombination3'));
            end
            
        end % end of checkDistributionLinkCombination.
        
        function sglme = checkForBadlyScaledPLWeights(sglme,sqrtdiagW)
% Takes a N-by-1 vector of the square root of the weights from PL
% iterations. We check if the weights have an enormous range using a
% technique similar to glmfit. If we detect bad scaling, we issue a warning
% message. To avoid repeated warning messages like this, we issue the
% warning message only if sglme.HaveWarnedAboutBadlyScaledPLWeights is
% false. The returned object may have HaveWarnedAboutBadlyScaledPLWeights
% set to true if we just threw a warning.

            % If the weights have an enormous range, we won't be able to do
            % PL very well.  The prior weights may be bad, or the fitted
            % mu's may have too wide a range, which is probably because the
            % data do as well, or because the link function is trying to go
            % outside the distribution's support. This technique is adapted
            % from glmfit.
            dataClass = class(sqrtdiagW);
            wtol = max(sqrtdiagW)*eps(dataClass)^(2/3);
            t = (sqrtdiagW < wtol);
            if any(t)
                t = t & (sqrtdiagW ~= 0);
                if any(t)   
                    warned = sglme.HaveWarnedAboutBadlyScaledPLWeights;
                    if ~warned
                        warning(message('stats:classreg:regr:lmeutils:StandardGeneralizedLinearMixedModel:BadlyScaledPLWeights'));
                        sglme.HaveWarnedAboutBadlyScaledPLWeights = true;
                    end
                end
            end

        end % end of checkForBadlyScaledPLWeights.        
    end
    
% Public, static, hidden validation methods.
    methods (Static,Access=public,Hidden=true)
        
        function linkStruct = validateLink(linkSpec,fitMethod)
%linkStruct = validateLink(sglme,linkSpec,fitMethod) takes a string, number
%   or structure linkSpec that specifies the link function g, validates it
%   and returns the validated link information as a structure linkStruct
%   which looks like this:
%
%        Name                 Name of the link function.
%        Link                 The function that defines g.
%        Derivative           Derivative of g.
%        SecondDerivative     Second derivative of g.
%        Inverse              Inverse of g.
%
%   We don't need to fill in the SecondDerivative info if fitMethod is
%   'MPL' or 'REMPL' but we need the second derivatives for 'Laplace',
%   'ApproximateLaplace' and 'Quadrature'.
%
%   If linkSpec is a structure, it must have the fields Link, Derivative
%   and Inverse. If fitMethod is *not* 'mpl' or 'rempl', it must also have
%   the field SecondDerivative.

            % (1) linkSpec must be a string, number or structure.
            isastring = internal.stats.isString(linkSpec);
            isanumber = isnumeric(linkSpec) & isreal(linkSpec) & isscalar(linkSpec);
            isastruct = isstruct(linkSpec);
            if ~(isanumber || isastruct || isastring)
                error(message('stats:classreg:regr:lmeutils:StandardGeneralizedLinearMixedModel:BadLinkSpec'));
            end            
            
            % (2) If a string, convert linkSpec to lower case.
            if isastring
                linkSpec = lower(linkSpec);
            end
            
            % (3) Ensure that fitMethod is a string.
            assert(internal.stats.isString(fitMethod));                        
            
            % (4) Helper function stattestlink gets g, g' and ginv:
            % [link,dlink,ilink] = stattestlink(linkArg,dataClass) 
            [linkStruct.Link,linkStruct.Derivative,linkStruct.Inverse] =...
                dfswitchyard('stattestlink',linkSpec,'double');

            % (5) Get the name of the link specified by linkSpec. Use the 
            % string 'custom' if linkSpec is a structure. We know for sure
            % that linkSpec is either a number, a structure or a string.
            if isanumber
                if linkSpec == 0 % equivalent to the log link
                    linkName = 'log';
                else
                    linkName = 'power';
                    linkExponent = linkSpec;
                end                
            elseif isastruct
                linkName = 'custom'; 
            elseif isastring
                % linkSpec is already lower case.
                linkName = linkSpec;
            end

            % (6) Add second derivative info depending on the linkSpec.            
            switch lower(linkName)
                case 'identity'
                    linkStruct.SecondDerivative = @(mu) zeros(size(mu));
                    linkStruct.Name = 'identity';
                case 'log'
                    linkStruct.SecondDerivative = @(mu) -1./(mu.^2);
                    linkStruct.Name = 'log';
                case 'logit'
                    linkStruct.SecondDerivative = @(mu) (2*mu-1)./((mu.*(1-mu)).^2);
                    linkStruct.Name = 'logit';
                case 'probit'
                    linkStruct.SecondDerivative = @(mu) norminv(mu)./((normpdf(norminv(mu))).^2);
                    linkStruct.Name = 'probit';
                case 'comploglog'
                    linkStruct.SecondDerivative = @(mu) -(1 + log(1-mu))./(((1-mu).*log(1-mu)).^2);
                    linkStruct.Name = 'comploglog';
                case {'loglog', 'logloglink'}
                    linkStruct.SecondDerivative = @(mu) -(1 + log(mu))./((mu.*log(mu)).^2);
                    linkStruct.Name = 'loglog';
                case 'reciprocal'
                    linkStruct.SecondDerivative = @(mu) 2./(mu.^3);
                    linkStruct.Name = 'reciprocal';
                case 'power' % linkExponent==0 (equivalent to 'log') has been weeded out already
                    linkStruct.SecondDerivative = @(mu) linkExponent*(linkExponent-1)*(mu.^(linkExponent-2));
                    linkStruct.Name = 'power';
                case 'custom'
                    % linkSpec is a struct. If fitmethod is 'mpl' or
                    % 'rempl', we don't need the second derivative info.                    
                    if ~any(strcmpi(fitMethod,{'mpl','rempl'}))
                        % Validate linkSpec.SecondDerivative and save it
                        % into linkStruct.SecondDerivative.
                        if ~isfield(linkSpec,'SecondDerivative')
                            error(message('stats:classreg:regr:lmeutils:StandardGeneralizedLinearMixedModel:MustSupplySecondDerivativeForML'));    
                        end
                        d2link = linkSpec.SecondDerivative;
                        if ischar(d2link) && ~isempty(which(d2link))
                            name = d2link; d2link = @(mu) feval(name,mu);
                        elseif ~isa(d2link,'function_handle') && ~isa(d2link,'inline')
                            error(message('stats:classreg:regr:lmeutils:StandardGeneralizedLinearMixedModel:BadLinkSpecSecondDerivative'));
                        end
                        linkStruct.SecondDerivative = d2link;                              
                    end
                    
                    % Copy 'Name' from linkSpec if you can.
                    if isfield(linkSpec,'Name')
                        linkStruct.Name = linkSpec.Name;
                    else
                        linkStruct.Name = 'custom';
                    end
                    
                otherwise
                    error(message('stats:stattestlink:UnrecognizedLink'));
            end  
            
        end % end of validateLink.
        
    end
    
% Protected methods related to fitting.
    methods (Access=protected)        
        
        function weights = getEffectiveObservationWeights(sglme)
%weights = getEffectiveObservationWeights(sglme) combines the prior weights 
%   with binomial size if required and gets the effective observation 
%   weights. For most distributions the prior weights are the same as
%   effective weights. The one exception is the binomial distribution for
%   which the prior weights need to be multiplied by the binomial size.
            
            if strcmpi(sglme.Distribution,'binomial') == 1
                weights = sglme.PriorWeights .* sglme.BinomialSize;                
            else
                weights = sglme.PriorWeights;                
            end
            
        end % end of getEffectiveObservationWeights.                
        
        function mu = initializeMuForPL(sglme)
%mu = initializeMuForPL(sglme) initializes the conditional mean of y given 
%   b for PL iterations. There are three possible situations:
%
%   (1) If sglme is a previously fitted model, we compute mu based on 
%   values of betaHat and bHat stored in the object.
%
%   (2) If sglme is an unfitted model and sglme.MuStart is empty, we know 
%   that sglme.y is already a valid value for the specified distribution.
%   So we can adjust sglme.y further to avoid boundary values. This mimics
%   the technique used in glmfit.
%
%   (3) If sglme is an unfitted model and sglme.MuStart is given, we simply
%   set mu to be equal to MuStart. This is a user supplied initial value
%   for mu.

            if (sglme.isFitToData == true)
                % Initialization using previously fitted values.
                X = sglme.X;
                Z = sglme.Z;
                delta = sglme.Offset;
                betaHat = sglme.betaHat;
                bHat = sglme.bHat;                
                etaHat = X*betaHat + Z*bHat + delta;
                mu = sglme.Link.Inverse(etaHat);                
            else
                if isempty(sglme.MuStart)
                    % Default initialization.
                    y = sglme.y;
                    N = sglme.BinomialSize;
                    switch lower(sglme.Distribution)
                        case 'poisson'
                            mu = y + 0.25;
                        case 'binomial'
                            mu = (N .* y + 0.5) ./ (N + 1);
                        case {'gamma' 'inverse gaussian'}
                            mu = max(y, eps(class(y))); % somewhat arbitrary
                        otherwise
                            mu = y;
                    end
                else
                    % User supplied value.
                    mu = sglme.MuStart;
                end
            end
            
            % Ensure that mu respects distribution specific bounds.
            mu = constrainMu(sglme,mu,sglme.Distribution);
            
            % All elements of mu must be in the domain of the link
            % function. This may not happen with non-standard choices of
            % distribution/link combinations. If any element of mu is not
            % in the domain of the link function, warn about this.
            tf = isMuStartFeasibleForLink(sglme,mu);
            if ~tf
                warning(message('stats:classreg:regr:lmeutils:StandardGeneralizedLinearMixedModel:MuStartNotFeasible'));    
            end
            
            % Modify mu to ensure that all its elements are in the domain
            % of the link.
            mu = constrainMuForLink(sglme,mu,sglme.Link.Name);                        
            
        end % end of initializeMuForPL.            
        
        function tf = isMuStartFeasibleForLink(sglme,mu)
%tf = isMuStartFeasibleForLink(sglme,mu) returns true if the supplied or 
%   default value of 'MuStart' is sensible for the link and false
%   otherwise. The idea is to apply the link function to mu and see if we
%   get all real and finite values.

            % 1. Get the link.
            g = sglme.Link.Link;
            
            % 2. Get the linear predictor eta for mu.
            eta = g(mu);
            
            % 3. All elements of eta must be real and finite values.            
            if isreal(eta) && all(isfinite(eta))
                tf = true;
            else
                tf = false;
            end
            
        end % end of isMuStartFeasibleForLink.
        
        function [tfeta,tfmu] = isLinearPredictorFeasible(sglme,eta)
%[tfeta,tfmu] = isLinearPredictorFeasible(sglme,eta) checks the supplied
%   value of linear predictor eta after a PL or ML fit. The idea is to
%   create a reconstruction of eta by transforming it to mu, constraining
%   mu based on the distribution and then transforming back to eta. Output
%   tfeta is true if the reconstruction error in eta is small. Similarly,
%   output tfmu is true if the reconstruction error in mu is small.

            % 1. Get the link and its inverse.
            g    = sglme.Link.Link;
            ginv = sglme.Link.Inverse;
            
            % 2. Transform eta to mu. Impose distribution specific
            % constraints on mu to get muCon.
            mu = ginv(eta);
            muCon = constrainMu(sglme,mu,sglme.Distribution);
            
            % 3. Transform muCon to etaCon - a reconstruction ot eta.
            etaCon = g(muCon);
            
            % 4. Ideally, etaCon should be equal to eta and muCon should be
            % equal to mu. Measure errors in reconstructions of eta and mu
            % on a relative scale.
            errEta = max(abs(etaCon - eta))/max(max(abs(eta)),1);            
            errMu  = max(abs(muCon  - mu ))/max(max(abs(mu)),1);
            
            % 5. Use sqrt(eps) as a threshold to detect small error.
            tol = sqrt(eps);
            
            % 6. Check error in eta reconstruction.
            if ( errEta <= tol )
                % eta is feasible.
                tfeta = true;
            else
                tfeta = false;
            end
            
            % 7. Check error in mu reconstruction.
            if ( errMu <= tol )
                % mu is feasible.
                tfmu = true;
            else
                tfmu = false;
            end
            
        end % end of isLinearPredictorFeasible.
        
        function mu = constrainMuForLink(sglme,mu,linkname)
%mu = constrainMuForLink(sglme,mu,linkname) takes a potential value mu for
%   the conditional mean of y given b, a string linkname and returns a mu
%   that is in the domain of the link. Each link implicitly imposes 
%   constraints on mu as described below.
%
%   linkname can be one of 'identity','log','logit','probit','comploglog',
%   'loglog','reciprocal','power','custom'.
%
%   Link                    Possible values of mu
%   logit                   (0,1)
%   probit                  (0,1)
%   comploglog              (0,1)
%   loglog                  (0,1)
%   log                     (0,Inf)
%   power                   (0,Inf)
%   reciprocal              (0,Inf)
%   identity                (-Inf,Inf)

            % 1. Get TINY and BIG.
            TINY = sglme.MuBound.TINY;
            BIG  = sglme.MuBound.BIG; 
            
            % 2. Enforce link specific constraints.
            switch lower(linkname)              
                case {'logit','probit','comploglog','loglog'}                    
                    %isok = all(mu > 0 & mu < 1);
                    isok = all(mu > TINY & mu < 1 - TINY);
                    if ~isok
                        a = max(TINY,eps);
                        mu = sglme.constrainVector(mu,a,1-a);
                    end
                case {'log','power','reciprocal'}
                    %isok = all(mu > 0 & mu < Inf);
                    isok = all(mu > TINY & mu < BIG);
                    if ~isok
                        a = max(TINY,eps);
                        b = min(BIG,realmax);
                        mu = sglme.constrainVector(mu,a,b);
                    end                   
                case {'identity'}
                    %isok = all(mu > -Inf & mu < Inf);
                    isok = all(mu > -BIG & mu < BIG);
                    if ~isok                       
                        b = min(BIG,realmax);                        
                        mu = sglme.constrainVector(mu,-b,b);
                    end
            end

        end % end of constrainMuForLink.
        
        function mu = constrainMu(sglme,mu,distribution)
%mu = constrainMu(sglme,mu,distribution) takes a potential value mu for the
%   conditional mean of y given b, a string distribution and returns a mu
%   that respects the constraints on mu enforced by the specified
%   distribution. Each distribution implicitly enforces constraints on mu
%   as follows:
%
%   Distribution            Possible values of mu
%   Binomial                (0,1)
%   Poisson                 (0,Inf)
%   Gamma                   (0,Inf)
%   Inverse Gaussian        (0,Inf)
%   Normal                  (-Inf,Inf)

            % (1) Ensure that distribution is a string and mu is a numeric 
            % real column vector.            
            % assert(internal.stats.isString(distribution));
            % assert(iscolumn(mu) & isnumeric(mu));
            mu = real(mu);

           TINY = sglme.MuBound.TINY;
            BIG = sglme.MuBound.BIG; 
            
            % (2) Enforce distribution specific constraints.
            switch lower(distribution)              
                case 'binomial'                    
                    %isok = all(mu > 0 & mu < 1);
                    isok = all(mu > TINY & mu < 1 - TINY);
                    if ~isok
                        a = max(TINY,eps);
                        mu = sglme.constrainVector(mu,a,1-a);
                    end
                case {'poisson','gamma','inverse gaussian','inversegaussian'}
                    %isok = all(mu > 0 & mu < Inf);
                    isok = all(mu > TINY & mu < BIG);
                    if ~isok
                        a = max(TINY,eps);
                        b = min(BIG,realmax);
                        mu = sglme.constrainVector(mu,a,b);
                    end                   
                case {'normal','gaussian'}
                    %isok = all(mu > -Inf & mu < Inf);
                    isok = all(mu > -BIG & mu < BIG);
                    if ~isok                       
                        b = min(BIG,realmax);                        
                        mu = sglme.constrainVector(mu,-b,b);
                    end
            end

        end % end of constrainMu.
        
        function displayPLConvergenceInfo(sglme,iter,loglik,etaTilde,etaHat,diagW,kappa,showploptimizerdisplay,isPLconverged)            
%displayPLConvergenceInfo(sglme,iter,loglik,etaTilde,etaHat,diagW,kappa,showploptimizerdisplay,isPLconverged)
%   prints PL convergence info on screen. Description of inputs:
%
%   iter                   = PL iteration number
%   loglik                 = maximized log likelihood from the LME approximation using etaTilde
%   etaTilde               = previous estimate of linear predictor
%   etaHat                 = current estimate of linear predictor
%   diagW                  = a vector of iterative weights for PL based on etaTilde
%   kappa                  = convergence tolerance on the linear predictor for terminating PL 
%   showploptimizerdisplay = true if PL fitting shows optimizer display on screen and false otherwise
%   isPLconverged          = true if PL iterations have converged and false otherwise
%
%   For reference, here's what the convergence test looks like in 
%   fitUsingPL:
%
%   ( max(abs(etaHat - etaTilde)) <= max(kappa*max(abs(etaTilde)),kappa) )
%
%   Our definition of deltaEta is based on this test.

            % 1. What is the infinity norm of etaTilde.
            infnormEtaTilde = max(abs(etaTilde));
            
            % 2. What is the relative change in linear predictor?
            deltaEta = max(abs(etaHat - etaTilde))/max(max(abs(etaTilde)),1);
            
            % 3. What is the two norm of weights?
            twonormW = norm(diagW);
            
            % 4. If you take etaHat to muHat and back to etaHat, do you get
            % the same etaHat?
            muHat              = sglme.Link.Inverse(etaHat);
            muHat              = constrainMu(sglme,muHat,sglme.Distribution);
            etaHatRecon        = sglme.Link.Link(muHat);
            etaHatReconError   = max(abs(etaHat - etaHatRecon))/max(max(abs(etaHat)),1);
            
            % 5. Print header every 20 iterations and every iteration if PL
            % optimizer display is on.
            if ( showploptimizerdisplay || rem(iter,20) == 1 )  
                fprintf('\n');
                fprintf('  -----------------------------------------------------------------------------------\n');
                fprintf('  PL ITER      LOGLIK       ||ETA||    ||ERR: ETA||    ||W||    ||ERR: ETA->MU->ETA||\n');
                fprintf('  -----------------------------------------------------------------------------------\n');
            end
            
            % 6. Display convergence info.
            fprintf('%9d    %+6.3e    %06.3e    %06.3e    %06.3e       %06.3e\n',iter,loglik,infnormEtaTilde,deltaEta,twonormW,etaHatReconError);                        
            
            % 7. Also display final convergence info if PL has converged.
            % Here's an example:
            % PL converged: Relative change in linear predictor = 1.000e-03, PLTolerance = 1.000e-08
            if ( isPLconverged == true )                
                plconvergedstring    = getString(message('stats:classreg:regr:lmeutils:StandardGeneralizedLinearMixedModel:Display_PLIterationDetail1'));
                relchangestring      = getString(message('stats:classreg:regr:lmeutils:StandardGeneralizedLinearMixedModel:Display_PLIterationDetail2'));
                pltolstring          = getString(message('stats:classreg:regr:lmeutils:StandardGeneralizedLinearMixedModel:Display_PLIterationDetail3'));
                fprintf('\n');
                fprintf('%s: %s %6.3e, %s %6.3e\n',plconvergedstring,relchangestring,deltaEta,pltolstring,kappa);
            end
            
        end % end of displayPLConvergenceInfo.
                        
        function tf = showSummaryMessages(sglme)
%showSummaryMessages(sglme) returns true if summary messages need to be
%   printed to screen on top of iterative optimizer output. tf is true if
%   summary messages should be shown on screen and false otherwise. Here's
%   the logic:
%
%   1. If sglme.OptimizerOptions is empty then tf is false.
%
%   2. If sglme.OptimizerOptions is not empty then we query its 'Display'
%   field. If it contains 'off' or 'none' then tf is false. Otherwise tf is
%   true.

            % 1. Get optimizer options.
            optimizeroptions = sglme.OptimizerOptions;
            
            % 2. optimizeroptions is empty.
            if isempty(optimizeroptions)
                tf = false;
                return;
            end
            
            % 3. optimizeroptions is not empty. It may be a structure or an
            % object with a 'Display' field.
            if any(strcmpi(optimizeroptions.Display,{'off','none'}))
                tf = false;
            else
                tf = true;
            end
            
        end % end of showSummaryMessages.
        
        function checkFinalPLSolution(sglme,eta)
%checkFinalPLSolution(sglme,eta) is a utility function to ensure that the 
%   final PL solution is sensible.
            
            % 1. Check feasibility status of eta and the corresponding mu.
            [tfeta,tfmu] = isLinearPredictorFeasible(sglme,eta);
            
            % 2. Consider the following scenarios:
            %
            % (1) If both eta and mu are good, there is nothing to do. 
            %
            % (2) If both eta and mu are bad, we warn about this.
            %
            % (3) In some cases, elements of eta are unbounded but the 
            % corresponding mu is well defined (e.g., binomial distribution 
            % when some of the fitted values approach 0/1, Hauck-Donner 
            % effect etc). In such cases, tfeta is false but tfmu is true.
            % We throw a warning in such cases.
            %
            % (4) If tfeta is true but tfmu is false, that's a warning as
            % well.
            if     ( tfeta == true  && tfmu == true )
                % Good.
            elseif ( tfeta == false && tfmu == false )
                % Bad.
                warning(message('stats:classreg:regr:lmeutils:StandardGeneralizedLinearMixedModel:BadFinalPLSolution'));                
            else 
                % Warn about other 2 cases and attempt to continue.
                warning(message('stats:classreg:regr:lmeutils:StandardGeneralizedLinearMixedModel:BadFinalPLSolution'));
            end            

        end % end of checkFinalPLSolution.
        
        function [sglme,cause] = fitUsingPL(sglme,numIter,kappa)                        
%[sglme,cause] = fitUsingPL(sglme,numIter,kappa) attempts to fit the
%   specified GLME model via pseudo likelihood (PL). numIter is the maximum
%   number of PL iterations to perform and kappa is the relative
%   convergence tolerance on the linear predictor. sglme is the modified
%   StandardGeneralizedLinearMixedModel object with betaHat, bHat, sigmaHat,
%   thetaHat, loglikHat and Psi filled in. cause is an integer indicating
%   the reason for termination:
%
%   cause          Meaning
%     0            PL iterations terminated because of kappa tolerance.
%     1            PL iterations terminated because of iteration limit.
%
%   cause = 1 implies that PL iterations were not able to converge. If
%   sglme.FitMethod is 'MPL' or 'REMPL', we would like cause = 0. For other
%   values of sglme.FitMethod, PL may be used as an initialization method
%   and so cause = 1 is not necessarily bad.

            % 1.1 Get and validate numIter and kappa.
            if nargin < 2
                numIter = sglme.PLIterations;
                kappa = sglme.PLTolerance;
            elseif nargin < 3
                kappa = sglme.PLTolerance;
            end
            assert(isscalar(numIter) ...
                    & internal.stats.isIntegerVals(numIter,1));
            assert(isnumeric(kappa) & isreal(kappa) & isscalar(kappa));
            
            % 1.2 Set verbosity level based on sglme.OptimizerOptions. If
            % verbose is true, we show summary messages on top of the
            % optimizer display.
            verbose = showSummaryMessages(sglme);            
            
            % 2. Get X, y, Z, delta, effective obs weights, FitMethod and
            % distribution.
                       X = sglme.X;
                       y = sglme.y;
                       Z = sglme.Z;
                   delta = sglme.Offset;
                       w = getEffectiveObservationWeights(sglme);
               fitmethod = sglme.FitMethod;
            distribution = sglme.Distribution;
            
            % 3. Ensure that fitmethod is 'mpl' or 'rempl'. fitUsingPL may
            % be used to initialize ML based methods. So if fitmethod is
            % not 'mpl' or 'rempl', set it equal to 'mpl'.
            if ~any(strcmpi(fitmethod,{'mpl','rempl'}))
                fitmethod = 'mpl';
            end
            
            % 4. Get link function g, its inverse ginv, its derivative gp 
            % and variance function v.
               g = sglme.Link.Link;
            ginv = sglme.Link.Inverse;
              gp = sglme.Link.Derivative;
               v = sglme.VarianceFunction.VarianceFunction;
            
            % 5. Initialize muTilde - the conditional mean of y given b. 
            % Then get linear predictor etaTilde given muTilde.
             muTilde = initializeMuForPL(sglme);
            etaTilde = g(muTilde);
                        
            % 6. Do PL iterations.
            iter = 1;
            found = false;
            
            if ( verbose == true )
               fprintf('\n%s\n',getString(message('stats:classreg:regr:lmeutils:StandardGeneralizedLinearMixedModel:Display_StartingPL'))); 
            end
            
            while( found == false )    
                
                % 6.1 Compute muTilde from etaTilde and enforce
                % distribution specific constraints on muTilde.                
                muTilde = ginv(etaTilde);
                muTilde = constrainMu(sglme,muTilde,distribution);
                
                % 6.2 Compute pseudo data yp. (etaTilde - delta) enters the
                % definition of yp not etaTilde.
                yp = gp(muTilde).*(y - muTilde) + (etaTilde - delta);
                                               
                % 6.3 Fit the following LME model:
                % yp = X*beta + Z*b + epsilon
                % b ~ N(0,sigma^2*D)
                % cov(epsilon) = sigma^2 * inv(W).
                % W is a diagonal matrix.
                
                    % 6.3.1 Compute diagonal elements of W.
                    diagW = w./v(muTilde)./(gp(muTilde).^2);
                    
                    % 6.3.2 Get weighted data for StandardLinearMixedModel
                    % and check for badly scaled PL weights.
                    sqrtdiagW = sqrt(diagW);
                        sglme = checkForBadlyScaledPLWeights(sglme,sqrtdiagW);
                          ypw = sqrtdiagW.*yp;
                           Xw = bsxfun(@times,sqrtdiagW,X);
                           Zw = bsxfun(@times,sqrtdiagW,Z);
                    
                    % 6.3.3 There are two possible scenarios:
                    %
                    %   (a) For the first iteration, initialize an unfitted 
                    %   StandardLinearMixedModel object using specified
                    %   options for Optimizer, OptimizerOptions,
                    %   InitializationMethod and CheckHessian. If
                    %   sglme.DispersionFixed is true, also set
                    %   'ResidualStd' option to value 1.0 to fix the
                    %   residual standard deviation at 1.0 in the call to
                    %   StandardLinearMixedModel.
                    %
                    %   (b) For 2,3,4... iterations, just change the X, y, 
                    %   Z properties of previously fitted slme. Also set 
                    %   the InitializationMethod property to 'default' to 
                    %   use previously computed thetaHat as initial value 
                    %   in the current LME fit.                
                    if iter == 1                    
                        Psiw = sglme.Psi;
                        switch lower(fitmethod)
                            case 'mpl'
                                fitmethodw = 'ml';
                            case 'rempl'
                                fitmethodw = 'reml';
                        end
                        
                        % If ShowPLOptimizerDisplay is false, we don't want
                        % to display iterative PL optimizer display. We set
                        % the 'Display' field in optimizer options to
                        % 'off'. We do this only if optimizer options is
                        % not empty.
                        optimizeroptions = sglme.OptimizerOptions;
                        if ( sglme.ShowPLOptimizerDisplay == false )
                            if ~isempty(optimizeroptions)
                                optimizeroptions.Display = 'off';
                            end
                        end                        
                        
                        dofit = false;
                        dostats = false;
                        args = {'Optimizer',sglme.Optimizer,...
                            'OptimizerOptions',optimizeroptions,...
                            'InitializationMethod',sglme.InitializationMethod,...
                            'CheckHessian',sglme.CheckHessian};
                        if (sglme.DispersionFixed == true)
                            args = [args,{'ResidualStd',1.0}];  %#ok<AGROW>
                        end
                        slme = classreg.regr.lmeutils.StandardLinearMixedModel(Xw,...
                            ypw,Zw,Psiw,fitmethodw,dofit,dostats,args{:});
                    else
                        slme.X = Xw;
                        slme.y = ypw;
                        slme.Z = Zw;                        
                        slme.InitializationMethod = 'default';
                    end
                    
                    % 6.3.4 Refit slme using new data.
                    slme = refit(slme);
                
                % 6.4 Current betaHat, bHat.
                betaHat = slme.betaHat;
                   bHat = slme.bHat;
                
                % 6.5 Current etaHat.
                etaHat = X*betaHat + Z*bHat + delta;
                
                % 6.6 Check for convergence.                
                if ( max(abs(etaHat - etaTilde)) <= max(kappa*max(abs(etaTilde)),kappa) )
                    found = true;
                    cause = 0;
                elseif ( iter >= numIter )
                    found = true;
                    cause = 1;
                end
                
                % 6.7 Display convergence info on screen if requested.
                if ( verbose == true )
                    if ( found == true && cause == 0 )
                        isPLconverged = true;
                    else
                        isPLconverged = false;
                    end
                    displayPLConvergenceInfo(sglme,iter,slme.loglikHat,etaTilde,etaHat,diagW,kappa,sglme.ShowPLOptimizerDisplay,isPLconverged);                                        
                end
                
                % 6.8 Update etaTilde and increment iter.
                etaTilde = etaHat;
                    iter = iter + 1;
                
            end % end of while.
            
            % 7. Modify sglme object.
              sglme.betaHat  =  slme.betaHat;
              sglme.bHat     =  slme.bHat;
             sglme.DeltabHat =  slme.DeltabHat;
            sglme.sigmaHat   =  slme.sigmaHat;
              sglme.thetaHat =  slme.thetaHat;
              sglme.phiHat   =   [];
              sglme.Psi      =  slme.Psi;
              sglme.slme     =  slme;

            % 8. Log likelihood approximation for final pseudo data yp and 
            % original y. Note that slme.loglikHat contains the log
            % likelihood for ypw.
                                muHat = ginv(etaHat);
                                muHat = constrainMu(sglme,muHat,distribution);
                              gpmuHat = gp(muHat);
                                diagW = w./v(muHat)./(gpmuHat.^2);
            sglme.loglikHatPseudoData = slme.loglikHat ...
                                            + 0.5*sum(log(abs(diagW)));
                      sglme.loglikHat = sglme.loglikHatPseudoData ...
                                            + sum(log(abs(gpmuHat)));
          
            % 9. Ensure that final PL solution is good.
            checkFinalPLSolution(sglme,etaHat);
            
        end % end of fitUsingPL.
       
        function sglme = initstatsPL(sglme)
%initstatsPL - Makes a GLME object sglme that has been fitted using PL
%   ready for stats. We can assume that isFitToData is true and
%   sglme.FitMethod is either 'mpl' or 'rempl'. This would imply that
%   sglme.slme is non-empty.
            
            % 1. Ensure that sglme.slme is not empty.
            if isempty(sglme.slme)
                error(message('stats:classreg:regr:lmeutils:StandardGeneralizedLinearMixedModel:MustRefitFirst'));
            end
            
            % 2. Call initstats on sglme.slme.
            assert(sglme.slme.isFitToData);
            sglme.slme = initstats(sglme.slme);
            
            % 3. Fill in covbetaHat, covthetaHatlogsigmaHat and
            % covetaHatlogsigmaHat.
              sglme.covbetaHat           = sglme.slme.covbetaHat;
            sglme.covthetaHatlogsigmaHat = sglme.slme.covthetaHatlogsigmaHat;
            sglme.covetaHatlogsigmaHat   = sglme.slme.covetaHatlogsigmaHat;
            
            % 4. Fill in rankX.
            sglme.rankX = sglme.slme.rankX;
            
        end % end of initstatsPL.
        
        function sglme = computeAMDPreordering(sglme,theta0)
%sglme = computeAMDPreordering(sglme,theta0) computes the AMD preordering
%   of U'*C*U + I_q for later use in posterior mode estimation. Since the
%   pattern of non-zeros in U'*C*U + I_q is the same as the pattern of
%   non-zeros in U'*U + I_q, we simply work with U'*U + I_q for the
%   purposes of computing the preordering. theta0 is the value of theta
%   that is used to compute U.
            
            % 1. Get the U matrix.
            U = getU(sglme,theta0);

            % 2. Get I_q.
            q = sglme.q;
            Iq = spdiags(ones(q,1),0,q,q);
            
            % 3. Store preordering in the object.
            sglme.AMDOrder = amd(U'*U + Iq);            

        end % end of computeAMDPreordering.
        
        function tf = isModelPreInitialized(sglme)
%tf = isModelPreInitialized(sglme) takes a potentially pre-initialized
%   model sglme and checks if things that should be filled in after a fit
%   using PL or ML have been filled in or not. tf is true if the model is
%   OK and false otherwise. What do we check? The following items must be
%   filled in: betaHat, thetaHat, sigmaHat, bHat, DeltabHat.
            
            % 1. Check betaHat.
            betaHat = sglme.betaHat;
            okbetaHat = (length(betaHat) == sglme.p);
            
            % 2. Check thetaHat.
            thetaHat = sglme.thetaHat;
            okthetaHat = (length(thetaHat) == sglme.Psi.NumParametersExcludingSigma);
            
            % 3. Check sigmaHat.
            sigmaHat = sglme.sigmaHat;
            oksigmaHat = (length(sigmaHat) == 1);
            
            % 4. Check bHat.
            bHat = sglme.bHat;
            okbHat = (length(bHat) == sglme.q);
            
            % 5. Check DeltabHat.
            DeltabHat = sglme.DeltabHat;
            okDeltabHat = (length(DeltabHat) == sglme.q);
            
            % 6. Is all OK?
            tf = okbetaHat & okthetaHat & oksigmaHat & okbHat & okDeltabHat;

        end % end of isModelPreInitialized.
        
        function [sglme,cause,x0,fun,xHat] = fitUsingApproximateLaplace(sglme)
%[sglme,cause] = fitUsingApproximateLaplace(sglme) takes a model previously
%   initialized using PL and attempts to fit it using 'ApproximateLaplace'
%   method. cause is an integer code indicating the reason for termination.
%   See fitUsingML for a description of cause. The pre-initialized model
%   *must* have valid values of betaHat, thetaHat, sigmaHat, bHat and
%   DeltabHat. thetaHat and sigmaHat are used to initialize the beta
%   profiled log likelihood fitting whereas betaHat and DeltabHat are used
%   to initialize the joint posterior estimation of (beta,Deltab) for each
%   function evaluation.
%
%[sglme,cause,x0,fun,xHat] = fitUsingApproximateLaplace(sglme) also returns
%   the initial point x0, a function handle fun for the objective function
%   being minimized and the estimated solution xHat.

            % 1. Ensure that betaHat, thetaHat, sigmaHat, bHat and 
            % DeltabHat are not empty.
            assert(isModelPreInitialized(sglme) == true);            

            % 2. Get sigma0 and theta0 from the object.
            sigma0 = sglme.sigmaHat;
            theta0 = sglme.thetaHat;                        
            
            % 3. Get the length thetaHat.
            lenthetaHat = length(theta0);
                        
            % 4. Set these warning off and then back on.
            warnState = warning('query','all');
            warning('off','MATLAB:nearlySingularMatrix');
            warning('off','MATLAB:illConditionedMatrix');
            warning('off','MATLAB:singularMatrix');
            cleanupObj = onCleanup(@() warning(warnState));                                                              
            
            % 5. Here's a definition of objective function to be minimized
            % for ML based methods.
            %---------------------------------------------------------------------------------------------------------------------
            %                               Laplace                       |  ApproximateLaplace   |     Quadrature
            %      Distribution
            % 
            %    'binomial'/'poisson'      x = [beta;theta]                  x = [theta]              x = [beta;theta]          
            %                              m1(x)                              m2(x)                   m3(x)  
            %
            %         other                x = [beta;theta;log(sigma)]       x = [theta;log(sigma)]   x = [beta;theta;log(sigma)]
            %                              m1(x)                              m2(x)                   m3(x)
            %----------------------------------------------------------------------------------------------------------------------
            %
            % m1 = Laplacian negative log likelihood with or without sigma fixed.
            % m2 = Beta-profiled Laplacian negative log likelihood with or without sigma fixed.
            % m3 = Quadrature based negative log likelihood with or without sigma fixed.            
            % 5.1 Initialize parameters.
            if sglme.isSigmaFixed
                x0 = theta0;
            else
                x0 = [theta0;log(sigma0)];
            end
            % 5.2 Make objective function for minimization.
            fun = makeNegativeApproximateLaplacianLogLikelihood(sglme);
            % 5.3 Do minimization.
            verbose = showSummaryMessages(sglme);            
            if ( verbose == true )
                fprintf('\n%s\n',getString(message('stats:classreg:regr:lmeutils:StandardGeneralizedLinearMixedModel:Display_StartingApproximateLaplace')));    
            end
            [xHat,cause] = doMinimization(sglme,fun,x0);
            % 5.4 Extract solution from xHat.
            thetaHat = xHat(1:lenthetaHat);
            
            if isempty(thetaHat)
                thetaHat = zeros(0,1);
            end
            
            if sglme.isSigmaFixed
                sigmaHat = 1.0;
            else
                logsigmaHat = xHat(lenthetaHat + 1);
                sigmaHat = exp(logsigmaHat);
            end
            % 5.5 Compute loglikHat (including constant), betaHat,
            % DeltabHat and bHat.
            includeConst = true;
            [loglikHat,betaHat,DeltabHat,bHat] = loglikelihoodApproximateLaplace(sglme,thetaHat,sigmaHat,includeConst);
                                                                                           
            % 6. Store betaHat, bHat, DeltabHat, sigmaHat, thetaHat,
            % loglikHat, phiHat and Psi in the object.
              sglme.betaHat = betaHat;
                 sglme.bHat = bHat;
            sglme.DeltabHat = DeltabHat;
             sglme.sigmaHat = sigmaHat;
             sglme.thetaHat = thetaHat;
            sglme.loglikHat = loglikHat;
               sglme.phiHat = [];
                  sglme.Psi = setUnconstrainedParameters(sglme.Psi,sglme.thetaHat);
                  sglme.Psi = setSigma(sglme.Psi,sglme.sigmaHat);

        end % end of fitUsingApproximateLaplace.
        
        function [sglme,cause,x0,fun,xHat] = fitUsingLaplace(sglme)
%[sglme,cause] = fitUsingLaplace(sglme) takes a model previously
%   initialized using PL or 'ApproximateLaplace' and attempts to fit it
%   using 'Laplace'. cause is an integer code indicating the reason for
%   termination. See the method fitUsingML for more info on cause. The
%   pre-initialized model *must* have valid values of betaHat, thetaHat,
%   sigmaHat, bHat and DeltabHat. betaHat, thetaHat and sigmaHat are used
%   to initialize the log likelihood fitting whereas DeltabHat is used to
%   initialize the posterior estimation of Deltab for each function
%   evaluation.
%
%[sglme,cause,x0,fun,xHat] = fitUsingLaplace(sglme) also returns the
%   initial point x0, a function handle fun for the objective function
%   being minimized and the estimated solution xHat.

            % 1. Ensure that betaHat, thetaHat, sigmaHat, bHat and 
            % DeltabHat are not empty.
            assert(isModelPreInitialized(sglme) == true);            

            % 2. Get beta0, sigma0 and theta0 from the object.
             beta0 = sglme.betaHat;
            sigma0 = sglme.sigmaHat;
            theta0 = sglme.thetaHat;                        
            
            % 3. Get the length of betaHat and thetaHat.
             lenbetaHat = length(beta0);
            lenthetaHat = length(theta0);
                        
            % 4. Set these warning off and then back on.
            warnState = warning('query','all');
            warning('off','MATLAB:nearlySingularMatrix');
            warning('off','MATLAB:illConditionedMatrix');
            warning('off','MATLAB:singularMatrix');
            cleanupObj = onCleanup(@() warning(warnState));                        
            
            % 5. Here's a definition of objective function to be minimized
            % for ML based methods.
            %---------------------------------------------------------------------------------------------------------------------
            %                               Laplace                       |  ApproximateLaplace   |     Quadrature
            %      Distribution
            % 
            %    'binomial'/'poisson'      x = [beta;theta]                  x = [theta]              x = [beta;theta]          
            %                              m1(x)                              m2(x)                   m3(x)  
            %
            %         other                x = [beta;theta;log(sigma)]       x = [theta;log(sigma)]   x = [beta;theta;log(sigma)]
            %                              m1(x)                              m2(x)                   m3(x)
            %----------------------------------------------------------------------------------------------------------------------
            %
            % m1 = Laplacian negative log likelihood with or without sigma fixed.
            % m2 = Beta-profiled Laplacian negative log likelihood with or without sigma fixed.
            % m3 = Quadrature based negative log likelihood with or without sigma fixed.                            
            % 5.1 Initialize parameters.
            if sglme.isSigmaFixed
                x0 = [beta0;theta0];
            else
                x0 = [beta0;theta0;log(sigma0)];
            end
            % 5.2 Make objective function for minimization.
            fun = makeNegativeLaplacianLogLikelihood(sglme);
            % 5.3 Do minimization.
            verbose = showSummaryMessages(sglme);            
            if ( verbose == true )
                fprintf('\n%s\n',getString(message('stats:classreg:regr:lmeutils:StandardGeneralizedLinearMixedModel:Display_StartingLaplace')));    
            end
            [xHat,cause] = doMinimization(sglme,fun,x0);
            % 5.4 Extract solution from xHat.
             betaHat = xHat(1:lenbetaHat);            
            thetaHat = xHat(lenbetaHat+1 : lenbetaHat+lenthetaHat);
            
            if isempty(betaHat)
                betaHat = zeros(0,1);
            end
            if isempty(thetaHat)
                thetaHat = zeros(0,1);
            end
            
            if sglme.isSigmaFixed
                sigmaHat = 1.0;
            else
                logsigmaHat = xHat(lenbetaHat+lenthetaHat+1);
                sigmaHat = exp(logsigmaHat);
            end
            % 5.5 Compute loglikHat (including constant), DeltabHat, bHat.
            includeConst = true;
            [loglikHat,DeltabHat,bHat] = loglikelihoodLaplace(sglme,thetaHat,sigmaHat,betaHat,includeConst);
            
            % 6. Store betaHat, bHat, DeltabHat, sigmaHat, thetaHat,
            % loglikHat, phiHat and Psi in the object.
              sglme.betaHat = betaHat;
                 sglme.bHat = bHat;
            sglme.DeltabHat = DeltabHat;
             sglme.sigmaHat = sigmaHat;
             sglme.thetaHat = thetaHat;
            sglme.loglikHat = loglikHat;
               sglme.phiHat = [];
                  sglme.Psi = setUnconstrainedParameters(sglme.Psi,sglme.thetaHat);
                  sglme.Psi = setSigma(sglme.Psi,sglme.sigmaHat);            
                        
        end % end of fitUsingLaplace.
        
        function [sglme,cause] = fitUsingML(sglme)
%[sglme,cause] = fitUsingML(sglme) attempts to fit the specified GLME model 
%   via maximum likelihood (ML). sglme is the modified object with betaHat,
%   bHat, sigmaHat, thetaHat, phiHat, loglikHat and Psi filled in. cause is
%   an integer indicating the reason for termination. The interpretation is
%   the same as doMinimization method of StandardLinearLikeMixedModel.
%
%    cause           Meaning
%     0,1            ML estimation converged.
%      2             ML estimation did not converge.
            
            % 1. Prepare for ML by doing preliminary PL iterations. Do the
            % PL iterations only if sglme is not fit to data yet, otherwise
            % just extract beta0, sigma0 and theta0 from the object.
            if ( sglme.isFitToData == false )
                  numIter = sglme.InitPLIterations;
                    kappa = sglme.PLTolerance;
                [sglme,~] = fitUsingPL(sglme,numIter,kappa);
            end
            
            % 2. Attempt to precompute a AMD ordering if requested.
            if ( sglme.UseAMDPreordering == true )
                sglme = computeAMDPreordering(sglme,sglme.thetaHat);                
            end            
            
            % 3. Here's a definition of objective function to be minimized
            % for ML based methods.
            %---------------------------------------------------------------------------------------------------------------------
            %                               Laplace                       |  ApproximateLaplace   |     Quadrature
            %      Distribution
            % 
            %    'binomial'/'poisson'      x = [beta;theta]                  x = [theta]              x = [beta;theta]          
            %                              m1(x)                              m2(x)                   m3(x)  
            %
            %         other                x = [beta;theta;log(sigma)]       x = [theta;log(sigma)]   x = [beta;theta;log(sigma)]
            %                              m1(x)                              m2(x)                   m3(x)
            %----------------------------------------------------------------------------------------------------------------------
            %
            % m1 = Laplacian negative log likelihood with or without sigma fixed.
            % m2 = Beta-profiled Laplacian negative log likelihood with or without sigma fixed.
            % m3 = Quadrature based negative log likelihood with or without sigma fixed.            
            switch lower(sglme.FitMethod)                
                case 'laplace'
                    % 3.1 For sequential fitting, first fit using
                    % approximate Laplace and use that as initialization
                    % for Laplace.
                    if ( sglme.UseSequentialFitting == true )
                        [sglme,~] = fitUsingApproximateLaplace(sglme);
                    end
                    [sglme,cause,~,fun,xHat] = fitUsingLaplace(sglme);
                    
                case 'approximatelaplace'
                    [sglme,cause,~,fun,xHat] = fitUsingApproximateLaplace(sglme);
                    
                case 'quadrature'                    
                    % TODO
            end
            
            % 4. Do Hessian checks only if ML converged and only if asked.
            doHessianCheck = (cause == 0 || cause == 1) && (sglme.CheckHessian == true);
            if ( doHessianCheck == true )
                checkObjectiveFunctionHessianForML(sglme,fun,xHat);
            end
            
        end % end of fitUsingML.
        
        function checkObjectiveFunctionHessianForML(sglme,fun,xHat)
%checkObjectiveFunctionHessianForML(sglme,fun,xHat) takes a fitted 
%   Standard GLME object sglme, a function handle fun to the objective
%   function that was minimized in ML and the solution xHat and checks the
%   positive definiteness of the objective function Hessian at xHat.

            if ~isempty(xHat)                
                % 1. Get the Hessian of fun at xHat.
                wantRegularized = false;
                H = sglme.getHessian(fun,xHat,wantRegularized);
                
                % 2. Get the right message ID.
                switch lower(sglme.FitMethod)
                    case 'laplace'
                        msgID = 'stats:classreg:regr:lmeutils:StandardGeneralizedLinearMixedModel:Message_NotSPDHessian_Laplace';
                    case 'approximatelaplace'
                        msgID = 'stats:classreg:regr:lmeutils:StandardGeneralizedLinearMixedModel:Message_NotSPDHessian_ApproximateLaplace';
                    otherwise
                        error(message('stats:classreg:regr:lmeutils:StandardGeneralizedLinearMixedModel:BadMLFitMethod'));
                end
                
                % 3. Call the positive definiteness check utility function.
                sglme.checkPositiveDefinite(H,msgID);                
            end

        end % end of checkObjectiveFunctionHessianForML.                
        
        function sglme = initstatsML(sglme)
%   This method is responsible for the following tasks:
%
%   1. Fill in covbetaHat.
%   2. Fill in covthetaHatlogsigmaHat.
%   3. Fill in covetaHatlogsigmaHat.
%   4. Do positive definiteness checks on covthetaHatlogsigmaHat and 
%      covetaHatlogsigmaHat if CheckHessian is true.
%   5. This method *does not* fill in covbetaHatbHat the covariance needed
%      for joint inference on beta/b. This is because covbetaHatbHat can be
%      a large matrix and so it makes sense to access it via a method on an
%      as needed basis.

            % 1. First fill out the covariance matrices.
            switch lower(sglme.CovarianceMethod)                
                case 'conditional'                    
                    % 1.1 Approximate covariance of betaHat conditional on
                    % the estimated covariance parameters equal to the true
                    % values.
                    sglme.covbetaHat = approximateCovBetaHatAsFunctionOfParameters(sglme,...
                        sglme.betaHat,sglme.DeltabHat,sglme.thetaHat,sglme.sigmaHat);
                    
                    % 1.2 This should be filled out on an as-needed basis.
                    sglme.covetaHatlogsigmaHat = [];                    
                    
                    % 1.3 Don't need to fill this out for GLMEs.
                    sglme.covthetaHatlogsigmaHat = [];                    
                case 'jointhessian'
                    
                    % 1.1 Joint covariance of
                    % [betaHat;etaHat;log(sigmaHat)] using observed
                    % information matrix computed using the Laplace
                    % approximation.
                    [sglme.covbetaHat,sglme.covetaHatlogsigmaHat] = covBetaHatEtaHatLogSigmaHat(sglme);                    
                    
                    % 1.2 Don't need to fill this out for GLMEs.
                    sglme.covthetaHatlogsigmaHat = [];                       
            end
            
            % 2. Set the rank of X.
            sglme.rankX = rank(sglme.X);
            
        end
        
    end
          
% Protected methods to "make" the various likelihood functions.
    methods (Access=protected)
        
        function fun = makeNegativeApproximateLaplacianLogLikelihood(sglme)
%fun = makeNegativeApproximateLaplacianLogLikelihood(sglme) returns the 
%   negative approximate Laplacian log likelihood as a function of model
%   parameters collected in the vector x. If sglme.isSigmaFixed is false 
%   then x = [theta;log(sigma)] and if sglme.isSigmaFixed is true then 
%   x = theta.

            fun = @f5;
            function y5 = f5(x)
                
                % 1. Extract theta from x.
                theta = x(1:sglme.Psi.NumParametersExcludingSigma,1);
                                
                % 2. Extract sigma from x or set it to its fixed value.
                if sglme.isSigmaFixed
                    sigma = sglme.sigmaFixed;
                else
                    logsigma = x(end);
                    sigma = exp(logsigma);
                end
                
                % 3. Evaluate approximate Laplacian log likelihood at
                % (theta,sigma) without the constant and negate the result.
                includeConst = false;
                L = loglikelihoodApproximateLaplace(sglme,theta,sigma,includeConst);                
                y5 = -1*L;
                
                % 4. Don't let y5 be equal to -Inf.
                y5 = max(-realmax,y5);
                
            end


        end % end of makeNegativeApproximateLaplacianLogLikelihood.
        
        function fun = makeNegativeLaplacianLogLikelihood(sglme)
%fun = makeNegativeLaplacianLogLikelihood(sglme) returns the negative 
%   Laplacian log likelihood as a function of model parameters collected in
%   the vector x. If sglme.isSigmaFixed is true then x = [beta;theta] and
%   if isSigmaFixed is false then x = [beta;theta;log(sigma)].
            
            fun = @f1;
            function y1 = f1(x)
                % 1. Extract beta and theta from x.
                 beta = x(1:sglme.p,1);
                theta = x(sglme.p + 1 : sglme.p + sglme.Psi.NumParametersExcludingSigma,1);                
                
                % 2. Extract sigma from x or set it to its fixed value.
                if sglme.isSigmaFixed
                    sigma = sglme.sigmaFixed;
                else
                    logsigma = x(end);
                    sigma = exp(logsigma);
                end
                
                % 3. Evaluate loglikelihoodLaplace at (beta,theta,sigma) 
                % without the constant and negate the result.
                includeConst = false;
                switch lower(sglme.EBMethod)
                    case 'default'
                        L = loglikelihoodLaplace(sglme,theta,sigma,beta,includeConst);  
                    otherwise
                        L = loglikelihoodLaplace2(sglme,theta,sigma,beta,includeConst);
                end
                y1 = -1*L;    
                
                % 4. Don't let y1 be equal to -Inf.
                y1 = max(-realmax,y1);
            end
            
        end % end of makeNegativeLaplacianLogLikelihoodAsAFunctionOfBetaTheta.        
        
        function fun = makeNegativeLaplacianLogLikelihoodNaturalParameters(sglme)            
%fun = makeNegativeLaplacianLogLikelihoodNaturalParameters(sglme) returns 
%   the negative Laplacian log likelihood as a function of model parameters
%   collected in the vector x using the Natural parameterization. If
%   sglme.isSigmaFixed is true then x = [beta;eta] and if isSigmaFixed is
%   false then x = [beta;eta;log(sigma)].
            
            fun = @f3;
            function y3 = f3(x)
                % 1. Extract beta and eta from x.
                beta = x(1:sglme.p,1);
                 eta = x(sglme.p + 1 : sglme.p + sglme.Psi.NumParametersExcludingSigma,1);                 
                 
                % 2. Extract sigma from x or set it to its fixed value.
                if sglme.isSigmaFixed
                    sigma = sglme.sigmaFixed;
                else
                    logsigma = x(end);
                    sigma = exp(logsigma);
                end
                
                % 3. Set sigma and eta in sglme.Psi.
                Psi = sglme.Psi;
                Psi = setSigma(Psi,sigma);
                Psi = setNaturalParameters(Psi,eta);
                
                % 4. Get theta from Psi.
                theta = getUnconstrainedParameters(Psi);
                
                % 5. Evaluate loglikelihoodLaplace at (beta,theta,sigma)
                % without the constant and negate the result.
                includeConst = false;
                switch lower(sglme.EBMethod)
                    case 'default'
                        L = loglikelihoodLaplace(sglme,theta,sigma,beta,includeConst);  
                    otherwise
                        L = loglikelihoodLaplace2(sglme,theta,sigma,beta,includeConst);
                end
                y3 = -1*L;    
                
                % 6. Don't let y3 be equal to -Inf.
                y3 = max(-realmax,y3);
                
            end
            
        end % end of makeNegativeLaplacianLogLikelihoodNaturalParameters.
        
    end
    
% Public helper methods. Defined public for use by unit tests.
    methods (Access=public,Hidden=true)
        
        function diagW = getDiagW(sglme,mu,w) % Get diagonal W matrix.
%diagW = getDiagW(sglme,mu,w) takes a M-by-1 vector mu and a M-by-1 vector 
%   w and gets the diagonal elements of the W matrix associated with the 
%   GLME model.
            
            % 1. Get derivative of link gp and variance function v.              
            gp = sglme.Link.Derivative;
             v = sglme.VarianceFunction.VarianceFunction;            
        
            % 2. Get diagonal elements of W.
            diagW = w./v(mu)./(gp(mu).^2);               

        end % end of getDiagW.

        function diagC = getDiagC(sglme,mu,w) % Get diagonal C matrix.
%diagC = getDiagC(sglme,mu,w) takes a M-by-1 vector mu and a M-by-1 vector 
%   w and gets the diagonal elements of the C matrix associated with the 
%   GLME model.

            % 1. Get derivative of link gp, 2nd derivative of link g2p, 
            % variance function v, derivative of variance function vp.
             gp = sglme.Link.Derivative;
            g2p = sglme.Link.SecondDerivative;            
              v = sglme.VarianceFunction.VarianceFunction;
             vp = sglme.VarianceFunction.Derivative;            
            
            % 2. Form the vector xi.
            gpmu = gp(mu);
             vmu = v(mu);            
              xi = (g2p(mu).*vmu)./gpmu + vp(mu);
            
            % 3. Get diagW.
            diagW = w./vmu./(gpmu.^2); 
            
            % 4. Get diagC.
                y = sglme.y;
            diagC = (((y-mu).*xi)./vmu + 1).*diagW;             

        end % end of getDiagC.
        
        function cloglik = conditionalLogLikelihood(sglme,mu,sigma,includeConst) % Conditional GLME log likelihood.
%cloglik = conditionalLogLikelihood(sglme,mu,sigma,includeConst) takes a
%   StandardGeneralizedLinearMixedModel object sglme, a N-by-1 mean vector
%   mu, a scalar sigma representing the square root of dispersion parameter
%   and returns the conditional log likelihood of the observed data given
%   mu and sigma. If includeConst is true, cloglik is the full conditional
%   log likelihood. If includeConst is false, terms in the full conditional
%   log likelihood that do not depend on mu or sigma may be dropped and a
%   partial conditional log likelihood may be returned. For optimization
%   purposes, includeConst should be false but for reporting purposes
%   includeConst should be true.
            
            % 1. Validate mu, sigma and includeConst.
            assert(iscolumn(mu) & size(mu,1) == sglme.N);
            assert(isscalar(sigma) & all(sigma > 0));
            assert(isscalar(includeConst) & islogical(includeConst));
            
            % 2. Get the current distribution.
            distribution = sglme.Distribution;                     

            % 3. Get effective observation weights and response vector y.
            w = getEffectiveObservationWeights(sglme);
            y = sglme.y;
            
            % 4. Conditional log likelihood log(P(y | mu, sigma, w)).
            switch lower(distribution)
                case {'normal','gaussian'}
                    sigma2 = sigma^2;
                    cloglik = -(0.5/sigma2)*sum(w.*((y-mu).^2)) - 0.5*sum(log((2*pi*sigma2)./w));
                case 'binomial'
                    cloglik = sum(w.*(y.*log(mu) + (1-y).*log(1-mu)));
                    if (includeConst == true)
                        const = sum(gammaln(w+1) - gammaln(w.*y + 1) - gammaln(w.*(1-y) + 1));
                        cloglik = cloglik + const;
                    end
                case 'poisson'
                    cloglik = sum(w.*(y.*log(mu) - mu)) + sum(w.*y.*log(w));
                    if (includeConst == true)
                        const = -sum(gammaln(w.*y + 1));
                        cloglik = cloglik + const;
                    end
                case 'gamma'
                    sigma2 = sigma^2;
                    cloglik = sum((log((y.*w)./(mu*sigma2)) - (y./mu)).*(w/sigma2)) ...
                        - sum(log(y)) - sum(gammaln(w/sigma2));
                case {'inverse gaussian','inversegaussian'}
                    sigma2 = sigma^2;
                    cloglik = -(0.5/sigma2)*sum((w.*((y-mu).^2))./((mu.^2).*y)) ...
                        + 0.5*sum(log(w./(2*pi*sigma2*(y.^3))));
            end
            
        end % end of conditionalLogLikelihood.        

    end
    
% Protected helper methods.
    methods (Access=protected)
        
        function diagF = getDiagF(sglme,mu) % Get diagonal F matrix.
%diagF = getDiagF(sglme,mu) takes a M-by-1 vector mu and gets the diagonal
%   elements of the F matrix associated with the GLME model.
            
            % 1. Get derivative of link gp.              
            gp = sglme.Link.Derivative;
            
            % 2. Get diagonal elements of F.
            diagF = gp(mu);               

        end % end of getDiagF.
  
        function [U,Lambda] = getU(sglme,theta) % Get U and Lambda.
%U = getU(sglme,theta) takes a value of unconstrained parameters theta and
%   returns the U matrix associated with this theta.
%
%[U,Lambda] = getU(sglme,theta) also returns the lower triangular Cholesky
%   factor of D matrix in sglme.Psi.
    
            % 1. Get N-by-q matrix Z.
            Z = sglme.Z;
            
            % 2. Get sglme.Psi, set current theta and set sigma = 1.
            Psi = sglme.Psi;
            Psi = setUnconstrainedParameters(Psi,theta);
            Psi = setSigma(Psi,1);
            
            % 3. Get lower triangular Cholesky factor of D matrix in Psi.
            % Lambda has size q-by-q where q = size(Z,2).
            Lambda = getLowerTriangularCholeskyFactor(Psi);
            
            % 4. Compute the U matrix. U is N-by-q.
            U = Z*Lambda;                
            
        end % end of getU.

        function M = getUtCUPlusIdentity(sglme,U,diagC)
%M = getUtCUPlusIdentity(sglme,U,diagC) takes a N-by-q matrix U, a N-by-1 
%   vector diagC containing the diagonal elements of the C matrix and
%   returns M = U'*C*U + Iq.

%   Note: M computed using M = (U'*bsxfun(@times,diagC,U) + Iq) is *much
%   slower* compared to M computed using:
%   
%   q = sglme.q;
%   N = sglme.N;
%   diagC = spdiags(diagC,0,N,N);
%   Iq = spdiags(ones(q,1),0,q,q);
%   M = ((U'*(diagC*U)) + Iq);
            
            % 1. Get N and q.
            q = sglme.q;
            N = sglme.N;
            
            % 2. Make diagC into a sparse matrix.
            diagC = spdiags(diagC,0,N,N);
            
            % 3. Make q-by-q sparse identity matrix Iq.
            Iq = spdiags(ones(q,1),0,q,q);

            % 4. Make M = U'*C*U + Iq.
            M = ((U'*(diagC*U)) + Iq);

        end % end of getUtCUPlusIdentity.
        
        function beta0 = initializeBeta(sglme) % Initialize beta for ML.
%beta0 = initializeBeta(sglme) initializes beta0 for normalized posterior 
%   mode computation. The initial value is based on the final estimate from 
%   PL iterations. If PL estimate is not available then beta0 is set to a 
%   vector of all zeros.          

            if isempty(sglme.betaHat)
                beta0 = zeros(sglme.p,1);
            else
                beta0 = sglme.betaHat;
            end
            
        end % end of initializeBeta.

        function Deltab0 = initializeDeltab(sglme) % Initialize Deltab for ML.
%Deltab0 = initializeDeltab(sglme) initializes Deltab0 for normalized
%   posterior mode computation. The initial value is based on the final
%   estimate from initial PL iterations. If PL estimate is not available
%   then Deltab0 is set to a vector of all zeros.            

            if isempty(sglme.DeltabHat)
                Deltab0 = zeros(sglme.q,1);
            else
                Deltab0 = sglme.DeltabHat;
            end
            
        end % end of initializeDeltab.
        
    end
    
% Protected methods related to computing various likelihoods.
    methods (Access=protected)        
       
        function [logliklap,b] = loglikelihoodLaplaceAsFunctionOfParameters(sglme,theta,sigma,beta,Deltab,includeConst) % Laplacian log likelihood parameterized by beta, theta, sigma and Deltab.
%[logliklap,b] = loglikelihoodLaplaceAsFunctionOfParameters(sglme,theta,sigma,beta,Deltab,includeConst) 
%   takes given values of theta, sigma, beta and Deltab (normalized b) and 
%   computes the Laplacian log likelihood. For Laplace approximation,
%   Deltab is the posterior mode and for approximate Laplace [beta,Deltab]
%   is the joint posterior mode treating beta as random with a flat prior.
%   includeConst is true to include constants in the conditional log
%   likelihood of y given b or false to potentially exclude them. It is
%   sensible to set includeConst to false during optimization and true for
%   final reporting purposes.            

            % 1. Get X, q and offset delta from sglme.
                X = sglme.X;
                q = sglme.q;
            delta = sglme.Offset;
                        
            % 2. Get U and Lambda.
            [U,Lambda] = getU(sglme,theta);            
                        
            % 3. Get the linear predictor using beta and Deltab.
            eta = X*beta + U*Deltab + delta;
            
            % 4. Compute the mean vector.
            ginv = sglme.Link.Inverse;
              mu = ginv(eta);
              mu = constrainMu(sglme,mu,sglme.Distribution);
              
            % 5. get diagC.
                w = getEffectiveObservationWeights(sglme);
            diagC = getDiagC(sglme,mu,w);

            % 6. Get R and S such that S*R'*R*S' = U'*C*U + eye(q) where R
            % and S are q by q. S is a permutation matrix: S*S' = eye(q).
            % Iq = spdiags(ones(q,1),0,q,q);
            % [R,status,~] = chol(U'*bsxfun(@times,diagC,U) + Iq);            
            M = getUtCUPlusIdentity(sglme,U,diagC);
            [R,status,~] = chol(M);
            
            % 7. Make sure status is zero.
            if (status ~= 0)
                error(message('stats:classreg:regr:lmeutils:StandardGeneralizedLinearMixedModel:BadLaplaceLogLikelihood'));
            end

            % 8. Compute conditional log likelihood.
            cloglik = conditionalLogLikelihood(sglme,mu,sigma,includeConst);
            
            % 9. Compute the Laplacian log likelihood.
            logliklap = cloglik - 0.5*(Deltab'*Deltab)/(sigma^2) - sglme.logAbsDetTriangular(R);            

            % 10. Compute unnormalized posterior mode estimate if required.
            if nargout > 1
                b = Lambda*Deltab;
            end

        end % end of loglikelihoodLaplaceAsFunctionOfParameters.                                
        
        function [rb,eta,mu] = normalizedrb(sglme,y,X,beta,U,Deltab,delta,w) % rb, eta and mu for Laplace.
%[rb,eta,mu] = normalizedrb(sglme,y,X,beta,U,Deltab,delta,w) takes y, X, 
%   beta, U, Deltab, offset delta and effective observation weights vector
%   w and computes the q-by-1 vector of non linear equations corresponding
%   to the normalized posterior mode computation. The current eta and mu
%   are also returned.

            % 1. Get the mean vector mu.
            ginv = sglme.Link.Inverse;
            eta = X*beta + U*Deltab + delta;
            mu = ginv(eta);

            % 2. Enforce distribution specific constraint on mu.
            mu = constrainMu(sglme,mu,sglme.Distribution);
            
            % 3. Compute diagF and diagW.
            diagF = getDiagF(sglme,mu);
            diagW = getDiagW(sglme,mu,w);
            
            % 4. Compute rb.
            rb = U'*(diagW.*diagF.*(y-mu)) - Deltab;
            
        end % end of normalizedrb.
        
        function [Deltab,cause] = normalizedPosteriorModeOfB(sglme,theta,sigma,beta,numIter,kappa,stepTol) % Newton line search for normalized posterior mode of b.
%[Deltab,cause] = normalizedPosteriorModeOfB(sglme,theta,sigma,beta,numIter,kappa,stepTol)
%   computes the normalized posterior mode of the conditional distribution
%   of random effects given the response and current values of theta, sigma
%   and beta. A line search Newton method with backtracking line search is
%   used. This may fail if the Jacobian of the relevant system of non
%   linear equations becomes singular. numIter is the maximum number of
%   Newton iterations to perform and kappa is a relative convergence
%   tolerance on the linear predictor of the GLME. stepTol is the absolute
%   tolerance on the step length.
%
%   Deltab is the final estimated normalized posterior mode and cause has 
%   the following interpretation:
%
%   cause          Meaning
%     0            Posterior mode estimation converged (kappa tolerance).
%     1            Step size became smaller than stepTol.
%     2            Iteration limit reached without convergence.

            % 1. Get y, X, Z, q, Offset and effective observation weights.
                y = sglme.y;
                X = sglme.X;
                Z = sglme.Z;
                q = sglme.q;
            delta = sglme.Offset;
                w = getEffectiveObservationWeights(sglme);
            
            % 2. Get sglme.Psi, set current theta and set sigma = 1.
            Psi = sglme.Psi;
            Psi = setUnconstrainedParameters(Psi,theta);
            Psi = setSigma(Psi,1);
            
            % 3. Get lower triangular Cholesky factor of D matrix in Psi.
            % Lambda has size q-by-q where q = size(Z,2).
            Lambda = getLowerTriangularCholeskyFactor(Psi);
            
            % 4. Compute the U matrix. U is N-by-q.
            U = Z*Lambda;

            % 5. Compute initial value of Deltab.
            Deltab = initializeDeltab(sglme);
            
            % 6. Initial rb, eta and mu.
            [rb,eta,mu] = normalizedrb(sglme,y,X,beta,U,Deltab,delta,w);   
            
            % 7. Sufficient decrease constants. These are reasonable values
            % that should not need to be changed.
            rho = 0.5;
             c1 = 1e-5;
                         
            % 8. Do Newton iterations.
             iter = 0;
            found = false;             
            while( found == false )                                                              
                
                % 9.1 Get diagC.
                diagC = getDiagC(sglme,mu,w);

                % 9.2 Get Jb = -(U'*C*U + eye(q))
                %Iq = spdiags(ones(q,1),0,q,q);
                %Jb = -(U'*bsxfun(@times,diagC,U) + Iq);                
                Jb = -1*getUtCUPlusIdentity(sglme,U,diagC);                
                
                % 9.3 Solve for normalized search direction Deltap. Get R
                % and S such that S*R'*R*S' = -Jb where R and S are q by q
                % and S is a permutation matrix S*S' = eye(q).
                [R,status,S] = chol(-Jb);
                if (status ~= 0)
                    % Cholesky factorization didn't work.
                    Deltap = -(Jb \ rb);
                else
                    Deltap = S * ( R \ ( R' \ (S'*rb) ) );
                end
                
                % 9.4 Determine step length.
                     alpha = 1.0;
                foundalpha = false;                
                    % 9.4.1 Quantities that appear in sufficient decrease.
                    % We want:
                    % rbnew'*rbnew <= rb'*rb + 2*c1*alpha*(Deltap'*Jb'*rb)
                    % or
                    % rbnew'*rbnew <= term1  + 2*c1*alpha* term2
                    term1 = rb'*rb;
                    term2 = Deltap'*Jb'*rb;
                while ( foundalpha == false )                    
                    % 9.4.2 Compute new Deltab.
                    Deltabnew = Deltab + alpha*Deltap;
                    
                    % 9.4.3 Compute rbnew.
                    [rbnew,etanew,munew] = normalizedrb(sglme,y,X,beta,U,Deltabnew,delta,w);
                    
                    % 9.4.4 Check sufficient decrease.
                    % isok = (rbnew'*rbnew - (1-2*c1*alpha)*(rb'*rb)) <= 0;
                    isok = (rbnew'*rbnew - (term1  + 2*c1*alpha*term2)) <= 0;
                                        
                    % 9.4.5 Termination of step selection.
                    if ( isok == true )
                        foundalpha = true;
                    elseif ( alpha <= stepTol )
                        foundalpha = true;
                        warning(message('stats:classreg:regr:lmeutils:StandardGeneralizedLinearMixedModel:Message_PosteriorModeLineSearch'));
                    else
                        alpha = rho*alpha;
                    end           
                end
                                
                % 9.5 Check convergence on linear predictor.
                if ( max(abs(etanew - eta)) <= max(kappa*max(abs(eta)),kappa) )
                    found = true;
                    cause = 0;
                elseif ( norm(alpha*Deltap) <= stepTol )
                    found = true;
                    cause = 1;
                elseif ( iter >= numIter )
                    found = true;
                    cause = 2;
                end
                
                % 9.6 Update Deltab, rb, eta, mu and iter.
                Deltab = Deltabnew;
                    rb = rbnew;
                   eta = etanew;
                    mu = munew;                
                  iter = iter + 1;                
            end
            
        end % end of normalizedPosteriorModeOfB.
        
        function [logliklap,Deltab,b] = loglikelihoodLaplace(sglme,theta,sigma,beta,includeConst) % Laplacian log likelihood parameterized by beta, theta and sigma.
%[logliklap,Deltab,b] = loglikelihoodLaplace(sglme,theta,sigma,beta,includeConst) 
%   computes the Laplacian log likelihood of a GLME model as a function of
%   theta, sigma and beta. If includeConst is true, the full likelihood is
%   used in computing the conditional log likelihood. If includeConst is
%   false, constants in the conditional log likelihood not depending on
%   theta, sigma and beta may be omitted. For optimization, includeConst
%   should be false and for final reporting includeConst should be true.
%   Additional output arguments include the normalized posterior mode
%   Deltab and its unnormalized counterpart b.
            
            % 1. Get normalized posterior mode of random effects.            
            numIter = sglme.EBOptions.MaxIter;
              kappa = sglme.EBOptions.TolFun;
            stepTol = sglme.EBOptions.TolX;
            [Deltab,cause] = normalizedPosteriorModeOfB(sglme,theta,sigma,beta,numIter,kappa,stepTol);            
            if ( cause == 2 )
                warning(message('stats:classreg:regr:lmeutils:StandardGeneralizedLinearMixedModel:Message_PosteriorModeLineSearchNonConvergence'));
            end            
            
            % 2. Use Deltab to compute Laplacian log likelihood.
            [logliklap,b] = loglikelihoodLaplaceAsFunctionOfParameters(sglme,theta,sigma,beta,Deltab,includeConst);            
            
        end % end of loglikelihoodLaplace.               
                                                
        function [logliklap,Deltab,b] = loglikelihoodLaplace2(sglme,theta,sigma,beta,includeConst)
%[logliklap,Deltab,b] = loglikelihoodLaplace2(sglme,theta,sigma,beta,includeConst) 
%   computes the Laplacian log likelihood of a GLME model as a function of
%   theta, sigma and beta. If includeConst is true, the full likelihood is
%   used in computing the conditional log likelihood. If includeConst is
%   false, constants in the conditional log likelihood not depending on
%   theta, sigma and beta may be omitted. For optimization, includeConst
%   should be false and for final reporting includeConst should be true.
%   Additional output arguments include the normalized posterior mode
%   Deltab and its unnormalized counterpart b.
            
            % 1. Get normalized posterior mode of random effects.
            [Deltab,rDeltab,JDeltab,cause] = normalizedPosteriorModeOfB2(sglme,theta,sigma,beta);     
            if ( cause == 2 )
                warning(message('stats:classreg:regr:lmeutils:StandardGeneralizedLinearMixedModel:Message_PosteriorModeLineSearchNonConvergence'));
            end            
            
            % 2. JDeltab = U'*C*U + eye(q). Get R and S such that
            % S*R'*R*S' = U'*C*U + eye(q) where R and S are q by q. S is a
            % permutation matrix: S*S' = eye(q).
            if ~issparse(JDeltab)
                JDeltab = sparse(JDeltab);
            end
            [R,status,~] = chol(JDeltab);                        
            
            % 3. Make sure status is zero.
            if (status ~= 0)
                error(message('stats:classreg:regr:lmeutils:StandardGeneralizedLinearMixedModel:BadLaplaceLogLikelihood'));
            end

            % 4. Get X, Z and Offset from sglme.
                X = sglme.X;
                Z = sglme.Z;            
            delta = sglme.Offset;
            
            % 5. Get sglme.Psi, set current theta and set sigma = 1.
            Psi = sglme.Psi;
            Psi = setUnconstrainedParameters(Psi,theta);
            Psi = setSigma(Psi,1);
            
            % 6. Get lower triangular Cholesky factor of D matrix in Psi.
            % Lambda has size q by q where q = size(Z,2).
            Lambda = getLowerTriangularCholeskyFactor(Psi);
            
            % 7. Compute the U matrix. U is N by q.
            U = Z*Lambda;
                        
            % 8. Get the linear predictor using beta and Deltab.
            eta = X*beta + U*Deltab + delta;
            
            % 9. Compute the mean vector.
            ginv = sglme.Link.Inverse;
              mu = ginv(eta);
              mu = constrainMu(sglme,mu,sglme.Distribution);
              
            % 10. Compute conditional log likelihood.
            cloglik = conditionalLogLikelihood(sglme,mu,sigma,includeConst);
            
            % 11. Compute the Laplacian log likelihood.
            logliklap = cloglik - 0.5*(Deltab'*Deltab)/(sigma^2) - sglme.logAbsDetTriangular(R);            

            % 12. Compute unnormalized posterior mode estimate if required.
            if nargout > 2
                b = Lambda*Deltab;
            end                                   
            
        end % end of loglikelihoodLaplace2.               
                                
        function [loglikapproxlap,beta,Deltab,b] = loglikelihoodApproximateLaplace(sglme,theta,sigma,includeConst)
%[loglikapproxlap,beta,Deltab,b] = loglikelihoodApproximateLaplace(sglme,theta,sigma,includeConst)  
%   computes the approximate Laplacian log likelihood of a GLME model with
%   beta profiled out. If includeConst is true, the full likelihood is used
%   in computing the conditional log likelihood. If includeConst is false,
%   constants in conditional log likelihood not depending on theta, sigma
%   and beta may be omitted. For optimization, includeConst should be false
%   and for final reporting includeConst should be true. Additional output
%   arguments include the joint posterior modes beta and Deltab of P(beta,b
%   | y,theta,sigma) with a flat prior on beta and the unnormalized random
%   effects vector b corresponding to the normalized random effects vector
%   Deltab.
            
            % 1. Get the joint posterior mode of beta/Deltab.    
            [beta,Deltab,r,J,cause] = normalizedPosteriorModeOfBetaB(sglme,theta,sigma);
            if ( cause == 2 )
                warning(message('stats:classreg:regr:lmeutils:StandardGeneralizedLinearMixedModel:Message_JointPosteriorModeLineSearchNonConvergence'));
            end
            
            % 2. Evaluate Laplacian log likelihood.
            [loglikapproxlap,b] = loglikelihoodLaplaceAsFunctionOfParameters(sglme,theta,sigma,beta,Deltab,includeConst);

            
        end % end of loglikelihoodApproximateLaplace.
        
    end
    
% Public hidden methods related to posterior mode estimation of b or 
% [beta;b]. Define public for use by unit tests.
    methods (Access=public,Hidden=true)
        
        function rfun = makerfunForB(sglme,theta,sigma,beta) 
%rfun = makerfunForB(sglme,theta,sigma,beta) makes rfunForB as a function
%   of normalized random effects vector Deltab for given values of theta,
%   sigma and beta.

            % 1. Get y, X, Z, Offset and effective observation weights.
                y = sglme.y;
                X = sglme.X;
                Z = sglme.Z;                
            delta = sglme.Offset;
                w = getEffectiveObservationWeights(sglme);
            
            % 2. Get sglme.Psi, set current theta and set sigma = 1.
            Psi = sglme.Psi;
            Psi = setUnconstrainedParameters(Psi,theta);
            Psi = setSigma(Psi,1);
            
            % 3. Get lower triangular Cholesky factor of D matrix in Psi.
            % Lambda has size q by q where q = size(Z,2).
            Lambda = getLowerTriangularCholeskyFactor(Psi);
            
            % 4. Compute the U matrix. U is N by q.
            U = Z*Lambda;
            
            % 5. Get AMD related info.
                useamd = sglme.UseAMDPreordering;
                     s = sglme.AMDOrder;
            methodcode = sglme.NewtonStepMethodCode;
            
            
            % 5. rfun for posterior mode finding as a function of Deltab.
            rfun = @f2;
            function [r,J,pn] = f2(Deltab)
                %[r,J,pn] = rfunForB(sglme,Deltab,theta,sigma,beta,y,X,U,delta,w,useamd,s,methodcode)                
                switch nargout
                    case 1
                        r = rfunForB(sglme,Deltab,theta,sigma,beta,y,X,U,delta,w,useamd,s,methodcode);
                    case 2
                        [r,J] = rfunForB(sglme,Deltab,theta,sigma,beta,y,X,U,delta,w,useamd,s,methodcode);
                    case 3
                        [r,J,pn] = rfunForB(sglme,Deltab,theta,sigma,beta,y,X,U,delta,w,useamd,s,methodcode);
                end
            end            
            
        end % end of makerfunForB.

        function rfun = makerfunForBetaB(sglme,theta,sigma) 
%rfun = makerfunForBetaB(sglme,theta,sigma) makes rfunForBetaB as a 
%   function x = [beta;Deltab] where beta is the p-by-1 fixed effects
%   vector and Deltab is the q-by-1 normalized random effects vector for
%   given values of theta and sigma.

            % 1. Get y, X, Z, Offset and effective observation weights.
                y = sglme.y;
                X = sglme.X;
                Z = sglme.Z;                
            delta = sglme.Offset;
                w = getEffectiveObservationWeights(sglme);
            
            % 2. Get sglme.Psi, set current theta and set sigma = 1.
            Psi = sglme.Psi;
            Psi = setUnconstrainedParameters(Psi,theta);
            Psi = setSigma(Psi,1);
            
            % 3. Get lower triangular Cholesky factor of D matrix in Psi.
            % Lambda has size q by q where q = size(Z,2).
            Lambda = getLowerTriangularCholeskyFactor(Psi);
            
            % 4. Compute the U matrix. U is N by q.
            U = Z*Lambda;
            
            % 5. Get p, q, N.
            p = sglme.p;
            q = sglme.q;
            N = sglme.N;
            
            % 6. Get AMD related info.
                useamd = sglme.UseAMDPreordering;
                     s = sglme.AMDOrder;
            methodcode = sglme.NewtonStepMethodCode;                    
            
            % 7. rfun for posterior mode finding as a function of Deltab.
            rfun = @f4;
            function [r,J,pn] = f4(x)                
            % x = [beta;Deltab].                
                  beta = x(1:p,1);                
                Deltab = x(p+1:p+q,1);                
                
                switch nargout
                    case 1
                        r = rfunForBetaB(sglme,Deltab,theta,sigma,beta,y,X,U,delta,w,useamd,s,methodcode,N);
                    case 2
                        [r,J] = rfunForBetaB(sglme,Deltab,theta,sigma,beta,y,X,U,delta,w,useamd,s,methodcode,N);
                    case 3
                        [r,J,pn] = rfunForBetaB(sglme,Deltab,theta,sigma,beta,y,X,U,delta,w,useamd,s,methodcode,N);
                end
            end         
            
        end % end of makerfunForBetaB.
        
    end

% Protected methods related to posterior mode estimation of b or [beta;b].    
    methods (Access=protected)
        
        function [r,J,pn] = rfunForB(sglme,Deltab,theta,sigma,beta,y,X,U,delta,w,useamd,s,methodcode)
%[r,J,pn] = rfunForB(sglme,Deltab,theta,sigma,beta,y,X,U,delta,w,useamd,s,methodcode) 
%   makes the normalized set of non linear equations for normalized
%   posterior mode estimation. r is the q-by-1 vector of normalized
%   nonlinear equations evaluated at the specified value of normalized b
%   vector in Deltab and given values of theta, sigma and beta. J is the
%   Jacobian of r and pn is the solution to: (J'*J)*pn = -J'*r. J and pn
%   are computed on an as-needed basis.
% 
%   y is the GLME response, X is the fixed effects design matrix, U is the 
%   U(theta) matrix evaluated at current theta, delta is the offset vector 
%   and w is a vector of effective observation weights. useamd is true if
%   an AMD ordering should be used for computing pn and false otherwise. s
%   ia a q-by-1 preordering to use if useamd is true. If useamd is false, s
%   can be empty. methodcode is the technique to use for Newton step
%   computation. methodcode = 1 attempts to use 'Cholesky' first and if
%   that doesn't work falls back on \. methodcode = 2 always uses \.
%
%%
% Given $\theta$, $\sigma$ and $\beta$, the posterior mode of random
% effects is estimated as the solution to the following system of nonlinear
% equations:
%
% $r(b) = \Delta^{-T} \frac{\partial}{\partial b} \log{P(y,b | \beta, \theta, \sigma^2)} = 0$
%
% If we define the normalized random effects $b^*$ as $b^* = \Delta b$ then
% the system of nonlinear equations can also be written like this:
%
% $r^*(b^*) = \frac{\partial}{\partial b^*} \log{P(y,b | \beta, \theta, \sigma^2)} = 0$
%
% In terms of the GLME model parameters, this system of equations can be
% written like this (ignoring $\sigma^2$ in the denominator):
%
% $r^*(b^*) = \left\{ U(\theta)^T W(\mu(\beta,b)) F(\mu(\beta,b)) (y - \mu(\beta,b)) - b^*\right\}$
%
% We try to solve this normalized system of nonlinear equations. The
% Jacobian of $r^*(b^*)$ is given by:
%
% $J^*(b^*) = -\left\{ U(\theta)^T C(\mu(\beta,b)) U(\theta) + I_q \right\}$
%
% It would be more convenient if $J^*(b^*)$ were positive definite. To
% accomplish this we can multiply the set of equations $r^*(b^*)$ by -1.
% With this modification, the new set of equations and their Jacobian can
% be written like this:
%
% $r^*(b^*) = -\left\{ U(\theta)^T W(\mu(\beta,b)) F(\mu(\beta,b)) (y - \mu(\beta,b)) - b^*\right\}$
%
% $J^*(b^*) = \left\{ U(\theta)^T C(\mu(\beta,b)) U(\theta) + I_q \right\}$
%
% For simplicity of notation, let us denote this new $r^*(b^*)$ by $r$ and
% this new $J^*(b^*)$ by $J$. The Newton direction is given by $J* pn =
% -r$. We can solve this system via a Sparse Cholesky factorization of $J$.
% Suppose,
%
% $[L,status,S] = chol(J,'lower','matrix')$ then if $status = 0$ then $S$
% will be a permutation matrix with $S*S^T = I_q$ and $L$ will be a lower
% triangular matrix such that:
%
% $L*L^T = S^T * J * S$ or 
%
% $S*L*L^T*S^T = J$
%
% Thus $pn$ is given by:
%
% $S*L*L^T*S^T* pn = -r$ or
%
% $L*L^T*S^T* pn = -S^T*r$ or
%
% $L^T*S^T* pn = -(L \,\, bks \,\, (S^T*r))$ where $bks$ is "backslash".
% 
% $S^T* pn = -( L^T \,\, bks \,\, (L \,\, bks \,\, (S^T*r)) )$
%
% $pn = -( S* ( L^T \,\, bks \,\, (L \,\, bks \,\, (S^T*r)) ))$
%
% Instead of requesting $S$ as a matrix, we may request it as a vector like
% this:
%
% $[L,status,s] = chol(J,'lower','vector')$. The vector $s$ is related to
% the matrix $S$ as follows:
%
% $L*L^T = S^T * J * S = J(s,s)$ and for any compatible matrices $M$ and
% $T$:
%
% $S^T M = M(s,:)$ and 
%
% $T*S = T(:,s)$
%
% Thus, using $s$ the solution $pn$ can be calculated like this:
%
% $S^T* pn = -( L^T \,\, bks \,\, (L \,\, bks \,\, (S^T*r)) )$
%
% $pn(s) = -( L^T \,\, bks \,\, (L \,\, bks \,\, r(s)) )$
%
% For the above call, ensure that $pn$ is a column vector.
%
% If the permutation vector $s$ has been precomputed, we can apply it like
% this:
%
% $[L,status] = chol(J(s,s),'lower')$ followed by:
%
% $pn(s) = -( L^T \,\, bks \,\, (L \,\, bks \,\, r(s)) )$
%
%
% For example:
%
% [L,status,S] = chol(J,'lower','matrix');
%
% pn1 = -(S*(L'\(L\(S'*r))));
%
% [L,status,s] = chol(J,'lower','vector');
%
% pn2(s) = -(L'\(L\r(s)));
%
% [L,status] = chol(J(s,s),'lower');
%
% pn3(s) = -(L'\(L\r(s)));
%
% max(abs(pn1-pn2'))
% 
% ans =
% 
%      0
% 
% max(abs(pn1-pn3'))
% 
% ans =
% 
%      0

            % 1. Get the mean vector mu.
            ginv = sglme.Link.Inverse;
            eta = X*beta + U*Deltab + delta;
            mu = ginv(eta);

            % 2. Enforce distribution specific constraint on mu.
            mu = constrainMu(sglme,mu,sglme.Distribution);
            
            % 3. Compute diagF and diagW.
            diagF = getDiagF(sglme,mu);
            diagW = getDiagW(sglme,mu,w);
            
            % 4. Compute r.
            r = -(U'*(diagW.*diagF.*(y-mu)) - Deltab);
            
            % 5. Compute J if required.
            if nargout > 1            
                % 5.1 Compute diagC.
                diagC = getDiagC(sglme,mu,w);
                
                % 5.2 Compute J.                
                 J = getUtCUPlusIdentity(sglme,U,diagC);                 
            end

            % 6. Compute pn if required.
            if nargout > 2                                
                % 6.1 Solve J*pn = -r. Note that:
                %
                % J = U'*C*U + eye(q)                                
                %
                % is alrady positive definite.                
                pn = sglme.computeNewtonStepForB(J,r,useamd,s,methodcode);                
            end
            
        end % end of rfunForB.

        function [r,J,pn] = rfunForBetaB(sglme,Deltab,theta,sigma,beta,y,X,U,delta,w,useamd,s,methodcode,N)
%[r,J,pn] = rfunForBetaB(sglme,Deltab,theta,sigma,beta,y,X,U,delta,w) makes 
%   the normalized set of non linear equations for joint posterior mode
%   estimation of beta and b assuming a flat prior on beta. r is the
%   (p+q)-by-1 vector of normalized nonlinear equations evaluated at the
%   specified value of normalized b vector in Deltab, beta and given values
%   of theta and sigma. J is the Jacobian of r and pn is the solution to:
%   (J'*J)*pn = -J'*r. J and pn are computed on an as-needed basis.
% 
%   y is the GLME response, X is the fixed effects design matrix, U is the 
%   U(theta) matrix evaluated at current theta, delta is the offset vector 
%   and w is a vector of effective observation weights.
%
%%
% Given $\theta$, $\sigma$, the joint posterior mode of $(\beta,b)$ is 
% estimated as the solution to the following system of nonlinear equations:
%
% $r_{\beta}(\beta,b) = \frac{\partial}{\partial \beta} \log{P(y | \beta, b, \theta, \sigma^2)} = 0$
%
% $r_{b}(\beta,b) = \Delta^{-T} \frac{\partial}{\partial b} \log{P(y,b | \beta, \theta, \sigma^2)} = 0$
%
% If we define the normalized random effects $b^*$ as $b^* = \Delta b$ then
% the system of nonlinear equations can also be written like this:
%
% $r^*_{\beta}(\beta,b^*) = \frac{\partial}{\partial \beta} \log{P(y | \beta, b, \theta, \sigma^2)} = 0$
%
% $r^*_{b}(\beta, b^*) = \frac{\partial}{\partial b^*} \log{P(y,b | \beta, \theta, \sigma^2)} = 0$
%
% In terms of the GLME model parameters, this system of equations can be
% written like this (ignoring $\sigma^2$ in the denominator):
%
% $r^*_{\beta}(\beta,b^*) = X^T W(\mu(\beta,b)) F(\mu(\beta,b)) (y - \mu(\beta,b))$
%
% $r^*_{b}(\beta, b^*) = \left\{ U(\theta)^T W(\mu(\beta,b)) F(\mu(\beta,b)) (y - \mu(\beta,b)) - b^*\right\}$
%
% We try to solve this normalized system of nonlinear equations. The
% Jacobian of $[r^*_{\beta}(\beta,b^*);r^*_{b}(\beta, b^*)]$ is given by:
%
% $J^*(\beta,b^*) = -\left[ \begin{array}{cc} X^T C(\mu(\beta,b)) X & X^T C(\mu(\beta,b)) U \\ U^T C(\mu(\beta,b)) X & U^T C(\mu(\beta,b)) U + I_q \end{array}
% \right]$
%
% It is more convenient to make $J^*(\beta,b^*)$ positive definite. To that
% extent, we multiply the set of non linear equations by -1. The new system
% of non linear equations and their Jacobian is given by:
%
% $r^*_{\beta}(\beta,b^*) = -X^T W(\mu(\beta,b)) F(\mu(\beta,b)) (y - \mu(\beta,b))$
%
% $r^*_{b}(\beta, b^*) = -\left\{ U(\theta)^T W(\mu(\beta,b)) F(\mu(\beta,b)) (y - \mu(\beta,b)) - b^*\right\}$
%
% $J^*(\beta,b^*) = \left[ \begin{array}{cc} X^T C(\mu(\beta,b)) X & X^T C(\mu(\beta,b)) U \\ U^T C(\mu(\beta,b)) X & U^T C(\mu(\beta,b)) U + I_q \end{array}
% \right]$            
%
% Let us denote $r^*_{\beta}(\beta,b^*)$ by $r_{{\beta}}$, $r^*_{b}(\beta,
% b^*)$ by $r_{b}$ and $J^*(\beta,b^*)$ by $J$. Also let $r =
% [r_{\beta};r_b]$ and $J_b = U^T C(\mu(\beta,b)) U + I_q$.
%
% The "normalized" Newton step is given by (using new notation above): 
%
% $J * [p_{{\beta}};p_b] = -[r_{{\beta}};r_b]$ 
%
% Suppose,
%
% $[L,status,s] = chol(U^T C U + I_q,'lower','vector')$ with
%
% $S = I_q(:,s)$ 
%
% $Q1 = X^T  C  U  S  L^{-T}$
%
% $c_b = L^{-1}  S^T  r_b$ then
%
% $-p_{\beta} = (X^T C X - Q1 Q1^T)^{-1} (r_{{\beta}} - Q1 \, c_b)$
%
% $-p_{b} = S L^{-T} (c_b - Q1^T \, p_{\beta})$

            % 1. Get the mean vector mu.
            ginv = sglme.Link.Inverse;
            eta = X*beta + U*Deltab + delta;
            mu = ginv(eta);

            % 2. Enforce distribution specific constraint on mu.
            mu = constrainMu(sglme,mu,sglme.Distribution);
            
            % 3. Compute diagF and diagW.
            diagF = getDiagF(sglme,mu);
            diagW = getDiagW(sglme,mu,w);
            
            % 4. Compute r.
            p = sglme.p;
            q = sglme.q;            
            r = zeros(p+q,1);   
            
            rbeta = -(X'*(diagW.*diagF.*(y-mu)));
               rb = -(U'*(diagW.*diagF.*(y-mu)) - Deltab);            
                r(1:p) = rbeta;                        
            r(p+1:p+q) = rb;
            
            % 5. Compute J if required.
            if nargout > 1            
                % 5.1 Compute diagC.
                diagC = getDiagC(sglme,mu,w);
                
                % 5.2 Compute J.
                J = sparse(p+q,p+q);                                
                XtCX = X'*bsxfun(@times,diagC,X);
                %CU = bsxfun(@times,diagC,U);
                  CU = spdiags(diagC,0,N,N)*U;
                XtCU = X'*CU;
                UtCU = U'*CU;                
                  Iq = spdiags(ones(q,1),0,q,q);                
                UtCUIq = UtCU + Iq;  
                  
                      i1 = 1:p;
                      i2 = (p+1):(p+q);   
                J(i1,i1) = XtCX;
                J(i1,i2) = XtCU;
                J(i2,i1) = XtCU';
                J(i2,i2) = UtCUIq;                                             
            end

            % 6. Compute pn if required.
            if nargout > 2
                
                % 6.1 Solve J*pn = -r.                
                % pn = -(J \ r);                    
                pn = sglme.computeNewtonStepForBetaB(J,r,UtCUIq,XtCU,XtCX,rb,rbeta,p,q,Iq,useamd,s,methodcode);                                                              
                                    
            end
            
        end % end of rfunForBetaB.
        
        function [beta,Deltab,r,J,cause] = normalizedPosteriorModeOfBetaB(sglme,theta,sigma)
%[beta,Deltab,r,J,cause] = normalizedPosteriorModeOfBetaB(sglme,theta,sigma)  
%   computes the joint posterior mode of beta and Deltab for given values
%   of theta and sigma. beta and Deltab are the final posterior mode
%   estimates. r and J are the set of non linear equations corresponding to
%   joint posterior mode estimation and their Jacobian respectively
%   evaluated at final beta and Deltab. cause is an integer code indicating
%   the reason for termination:
%
%   cause          Meaning
%   0 or 1         Posterior mode estimation converged.
%     2            Iteration limit reached without convergence.

            % 1. Initialize beta and Deltab.
              beta0 = initializeBeta(sglme);
            Deltab0 = initializeDeltab(sglme);
            
            % 2. Set of non linear equations for beta/Deltab estimation.
            rfun = makerfunForBetaB(sglme,theta,sigma);
            
            % 3. Make initial point.
            x0 = [beta0;Deltab0];
            
            % 4. Call fsolve/nlesolve.
            if strcmpi(sglme.EBMethod,'fsolve')
                sglme.EBOptions.Jacobian = 'on';
                %[X,FVAL,EXITFLAG,OUTPUT,JACOB] = fsolve(FUN,X0,OPTIONS)
                [x,r,exitflag,~,J] = fsolve(rfun,x0,sglme.EBOptions);
                % 1  fsolve converged to a root.
                % 2  Change in X too small.
                % 3  Change in residual norm too small.
                % 4  Computed search direction too small.
                % 0  Too many function evaluations or iterations.
                % -1  Stopped by output/plot function.
                % -2  Converged to a point that is not a root.
                % -3  Trust region radius too small (Trust-region-dogleg).
                switch exitflag
                    case 1
                        cause = 0;
                    case {2,3,4}
                        cause = 1;
                    otherwise
                        cause = 2;
                end
            else
                [x,r,J,cause] = nlesolve(rfun,x0,...
                    'Options',sglme.EBOptions,'Method',sglme.EBMethod);                                  
            end
            
            % 5. Convert final x into beta, Deltab.
                 p = sglme.p;
                 q = sglme.q;            
              beta = x(1:p,1);             
            Deltab = x(p+1:p+q,1);                        
            
        end % end of normalizedPosteriorModeOfBetaB.

        function [Deltab,rDeltab,JDeltab,cause] = normalizedPosteriorModeOfB2(sglme,theta,sigma,beta)
%[Deltab,rDeltab,JDeltab,cause] = normalizedPosteriorModeOfB2(sglme,theta,sigma,beta)
%   computes the normalized posterior mode of the conditional distribution
%   of random effects given the response and current values of theta, sigma
%   and beta. A trust region based 2D subspace minimization method is used.
%   Deltab is the final estimated normalized posterior mode, rDeltab and
%   JDeltab are the normalized nonlinear equations and their Jacobian
%   evaluated at Deltab and cause is either 0 or 1 with the following
%   interpretation:
%
%   cause          Meaning
%   0 or 1         Posterior mode estimation converged.
%     2            Iteration limit reached without convergence.
            
            % 1. Initialize Deltab.
            Deltab0 = initializeDeltab(sglme);
                                    
            % 2. Set of nonlinear equations as a function of Deltab.
            rfun = makerfunForB(sglme,theta,sigma,beta);
            
            % 3. Call fsolve/nlesolve.
            if strcmpi(sglme.EBMethod,'fsolve')
                sglme.EBOptions.Jacobian = 'on';
                %[X,FVAL,EXITFLAG,OUTPUT,JACOB] = fsolve(FUN,X0,OPTIONS)
                [Deltab,rDeltab,exitflag,~,JDeltab] = fsolve(rfun,Deltab0,sglme.EBOptions);
                % 1  fsolve converged to a root.
                % 2  Change in X too small.
                % 3  Change in residual norm too small.
                % 4  Computed search direction too small.
                % 0  Too many function evaluations or iterations.
                % -1  Stopped by output/plot function.
                % -2  Converged to a point that is not a root.
                % -3  Trust region radius too small (Trust-region-dogleg).
                switch exitflag
                    case 1
                        cause = 0;
                    case {2,3,4}
                        cause = 1;
                    otherwise
                        cause = 2;
                end
            else
                [Deltab,rDeltab,JDeltab,cause] = nlesolve(rfun,Deltab0,...
                    'Options',sglme.EBOptions,'Method',sglme.EBMethod);            
            end
            
        end % end of normalizedPosteriorModeOfB2.
                        
    end
    
% Protected methods related to computing the covariance of betaHat, 
% [betaHat - beta;bHat - b] and [betaHat;etaHat;log(sigmaHat)]. 
    methods (Access=protected)
        
        function covbetaHat = approximateCovBetaHatAsFunctionOfParameters(sglme,beta,Deltab,theta,sigma)
%covbetaHat = approximateCovBetaHatAsFunctionOfParameters(sglme,beta,Deltab,theta,sigma) 
%   gets the covariance of betaHat as a function of beta, Deltab, theta and
%   sigma. Deltab is the normalized posterior mode of b calculated at beta,
%   theta and sigma. This covariance calculation is based on approximate
%   Hessian of the Laplacian log likelihood.

            % 1. Get X, p, q and Offset from sglme.
                X = sglme.X;
                p = sglme.p;
                q = sglme.q;
            delta = sglme.Offset;
            
            % 2. Get U matrix.
            U = getU(sglme,theta);                        
                        
            % 3. Get the linear predictor using beta and Deltab.
            eta = X*beta + U*Deltab + delta;
            
            % 4. Compute the mean vector.
            ginv = sglme.Link.Inverse;
              mu = ginv(eta);
              mu = constrainMu(sglme,mu,sglme.Distribution);
              
            % 5. Get diagC.
                w = getEffectiveObservationWeights(sglme);
            diagC = getDiagC(sglme,mu,w);

            
            % 6. Get C*U and C*X.
            CU = bsxfun(@times,diagC,U);
            CX = bsxfun(@times,diagC,X);
            
            % 7. Get X'*C*X and X'*C*U
            XtCX = X'*CX;
            XtCU = X'*CU;
            
            % 8. Get symmetrix matrix: M = U'*C*U + eye(q).
            Iq = spdiags(ones(q,1),0,q,q);
            M = U'*CU + Iq;
            M = 0.5*(M + M');
            
            % 9. Ensure that Cholesky factor of M exists, otherwise the
            % Laplace approximation is not well defined. Get R and S such
            % that S*R'*R*S' = M where R and S are q-by-q. S is a
            % permutation matrix: S*S' = eye(q).                                    
            [R,status,S] = chol(M);
                        
            % 10. Make sure status is zero.
            if (status ~= 0)
                error(message('stats:classreg:regr:lmeutils:StandardGeneralizedLinearMixedModel:BadLaplaceLogLikelihood'));
            end 
            
            % 11. Compute Q1.
            Q1 = (XtCU*S)/R;
            
            % 12. Get R1*R1'.
            R1R1t = XtCX - Q1*Q1';
            
            % 13. Get R1.
            [R1,status1] = chol(R1R1t,'lower');
            if (status1 ~= 0)
                % We know that R1R1t is positive definite. So inflate the
                % diagonal until Cholesky factorization succeeds.
                R1 = classreg.regr.lmeutils.covmats.CovarianceMatrix.singularLowerChol(R1R1t);
            end
            
            % 14. Get inverse of R1.
            R1inv = R1 \ eye(p);
            
            % 15. Get covariance of betaHat.
            covbetaHat = (sigma^2) * (R1inv'*R1inv);            
            
        end % end of approximateCovBetaHatAsFunctionOfParameters.
        
        function covbetaHatbHat = covBetaHatBHatAsFunctionOfOfParameters(sglme,beta,Deltab,theta,sigma)
%%
% The "joint" distribution of $\left[ \begin{array}{c} \beta \\ b \end{array}
% \right] \mid y, \hat{\theta},\hat{\sigma^2}$ is approximated by:
%
% $N\left( \left( \begin{array}{c} \hat{\beta} \\ \hat{b} \end{array}
% \right), \hat{\sigma^2} \, V(\hat{\theta},\hat{\sigma^2}) \right)$ where
% $V^{-1}$ can be written as:
%
% $V^{-1} = \left( \begin{array}{cc} X^T C(\mu(\hat{\beta},\hat{b})) X & X^T C(\mu(\hat{\beta},\hat{b})) Z \\ Z^T C(\mu(\hat{\beta},\hat{b})) X & Z^T C(\mu(\hat{\beta},\hat{b})) Z + D(\hat{\theta})^{-1} \end{array} \right)$
%
% If $U(\hat{\theta}) = Z \Lambda(\hat{\theta})$ with $D(\hat{\theta})^{-1} =
% \Delta(\hat{\theta})^T \Delta(\hat{\theta})$ and $\Lambda(\hat{\theta}) =
% \Delta(\hat{\theta})^{-1}$ then:
%
% $V^{-1} = \left( \begin{array}{cc} I_p & 0 \\ 0 & \Delta(\hat{\theta})^T \end{array} \right) \left( \begin{array}{cc} X^T C(\mu(\hat{\beta},\hat{b})) X & X^T C(\mu(\hat{\beta},\hat{b})) U(\hat{\theta}) \\ U(\hat{\theta})^T C(\mu(\hat{\beta},\hat{b})) X & U(\hat{\theta})^T C(\mu(\hat{\beta},\hat{b})) U(\hat{\theta}) + I_q \end{array} \right) \left( \begin{array}{cc} I_p & 0 \\ 0 & \Delta(\hat{\theta}) \end{array} \right)$
%
% $V^{-1} = K M K^T$
%
% $K = \left( \begin{array}{cc} I_p & 0 \\ 0 & \Delta(\hat{\theta})^T
% \end{array} \right)$
%
% $M = \left( \begin{array}{cc} X^T C(\mu(\hat{\beta},\hat{b})) X & X^T
% C(\mu(\hat{\beta},\hat{b})) U(\hat{\theta}) \\ U(\hat{\theta})^T C(\mu(\hat{\beta},\hat{b})) X & U(\hat{\theta})^T C(\mu(\hat{\beta},\hat{b})) U(\hat{\theta}) + I_q \end{array} \right)$
%
% If $M = L L^T$ then:
%
% $V^{-1} = K  L L^T K^T$
%
% $V = K^{-T} L^{-T} L^{-1} K^{-1}$
%
% $K^{-1} = G = \left( \begin{array}{cc} I_p & 0 \\ 0 & \Delta(\hat{\theta})^T
% \end{array} \right)^{-1} = \left( \begin{array}{cc} I_p & 0 \\ 0 & \Lambda(\hat{\theta})^T
% \end{array} \right)$
%
% $V = G^T L^{-T} L^{-1} G = T^T T$ where $T = L^{-1} G$.
%
% Output $covbetaHatbHat = \hat{\sigma^2} V$. For FitMethod equal to 'MPL' 
% or 'REMPL', we can replace $C(\mu(\hat{\beta},\hat{b}))$ by
% $W(\mu(\hat{\beta},\hat{b}))$ to get results consistent with LME.
            
            % 1. Get X, p, q, N and Offset from sglme.
                X = sglme.X;
                p = sglme.p;
                q = sglme.q;
                N = sglme.N;
            delta = sglme.Offset;
            
            % 2. Get U and Lambda matrices.
            [U,Lambda] = getU(sglme,theta);
                        
            % 3. Get the linear predictor using beta and Deltab.
            eta = X*beta + U*Deltab + delta;
            
            % 4. Compute the mean vector.
            ginv = sglme.Link.Inverse;
              mu = ginv(eta);
              mu = constrainMu(sglme,mu,sglme.Distribution);
              
            % 5. Get diagC and make it sparse diagonal matrix. For
            % FitMethod equal to 'mpl' or 'rempl', we replace the C matrix
            % by the W matrix to get the same results as from the final
            % fitted LME model from PL iterations.
            w = getEffectiveObservationWeights(sglme);
            switch lower(sglme.FitMethod)
                case {'mpl','rempl'}
                    diagC = getDiagW(sglme,mu,w); 
                otherwise
                    diagC = getDiagC(sglme,mu,w);
            end
            diagC = spdiags(diagC,0,N,N); 
            
            % 6. Get C*U and C*X.
            CU = diagC*U;
            CX = diagC*X;
            
            % 7. Get X'*C*X and X'*C*U
            XtCX = X'*CX;
            XtCU = X'*CU;
            
            % 8. Ensure that Cholesky factor of U'*C*U + eye(q) exists. If
            % not, the Laplace approximation is not well defined and hence
            % the joint covariance is also not well defined.
            Iq = spdiags(ones(q,1),0,q,q);         
            UtCUplusIq = U'*CU + Iq;
            [~,status,~] = chol(UtCUplusIq);
            if (status ~= 0)
                error(message('stats:classreg:regr:lmeutils:StandardGeneralizedLinearMixedModel:BadLaplaceLogLikelihood'));
            end
            
            % 9. Form the M matrix.
            M = sparse(p+q,p+q);
            i1 = 1:p;
            i2 = p+1:p+q;            
            M(i1,i1) = XtCX;
            M(i1,i2) = XtCU;
            M(i2,i1) = XtCU';
            M(i2,i2) = UtCUplusIq;
            
            % 10. Compute L such that M = L*L'.
            try 
                L = chol(M,'lower');
            catch ME %#ok<NASGU>
                L = classreg.regr.lmeutils.covmats.CovarianceMatrix.singularLowerChol(M);
            end
            
            % 11. Compute G = [I_p 0;0 Lambda^T].
            G = sparse(p+q,p+q);
            G(i1,i1) = eye(p);
            G(i2,i2) = Lambda';
            
            % 12. T = inv(L)*G.
            T = L \ G;
            
            % 13. Get sigma^2 * T'*T as the required covariance.
            covbetaHatbHat = (sigma^2)*(T'*T);
                                                                         
        end % end of covBetaHatBHatAsFunctionOfOfParameters.
        
        function C = covBetaHatBHat(sglme)
%covBetaHatBHat - Estimated covariance matrix of [betaHat-beta; bHat-b].          
%   C = covBetaHatBHat(sglme) returns the estimated covariance matrix of
%   [betaHat - beta; bHat - b].
 
            % Calculate the covariance if not already stored in object.                
            if isempty(sglme.covbetaHatbHat)
                  beta = sglme.betaHat;
                Deltab = sglme.DeltabHat;
                 theta = sglme.thetaHat;
                 sigma = sglme.sigmaHat;
                     C = covBetaHatBHatAsFunctionOfOfParameters(sglme,beta,Deltab,theta,sigma);
            else
                     C = sglme.covbetaHatbHat;
            end
            
        end % end of covBetaHatBHat.
        
        function [covbetaHat,covetaHatlogsigmaHat] = covBetaHatetaHatlogsigmaHatAsFunctionOfBetaThetaSigma(sglme,beta,theta,sigma)
%[covbetaHat,covetaHatlogsigmaHat] = covBetaHatetaHatlogsigmaHatAsFunctionOfBetaThetaSigma(sglme,beta,theta,sigma) 
%   computes the joint covariance of [betaHat;etaHat;log(sigmaHat)] using
%   the joint Hessian of the negative Laplacian log likelihood with respect
%   to [beta;eta;log(sigma)] evaluated at the [beta;eta;log(sigma)] values
%   corresponding to the supplied beta, theta and sigma.
            
            % 1. Make negative Laplacian log likelihood as a function of
            % [beta;eta] or [beta;eta;log(sigma)].
            fun = makeNegativeLaplacianLogLikelihoodNaturalParameters(sglme);
            
            % 2. Convert (theta,sigma) into (eta,sigma).
            Psi = sglme.Psi;
            Psi = setUnconstrainedParameters(Psi,theta);
            Psi = setSigma(Psi,sigma);            
            eta = getNaturalParameters(Psi);
            
            % 3. Does fun accept [beta;eta] or [beta;eta;log(sigma)]?
            if sglme.isSigmaFixed            
                % fun accepts [beta;eta].
                x = [beta;eta];                
            else                
                % fun accepts [beta;eta;log(sigma)].
                x = [beta;eta;log(sigma)];
            end
            
            % 4. Get the Hessian of fun w.r.t [beta;eta] or w.r.t
            % [beta;eta;log(sigma)] evaluated at x.
            try
                H = sglme.getHessian(fun,x);
            catch ME %#ok<NASGU>
                H = NaN(length(x));
            end
            
            % 5. If H is the Hessian of fun w.r.t [beta;eta] then append a
            % column and row of all zeros.
            if sglme.isSigmaFixed    
                n = size(H,1);
                H(n+1,:) = 0;
                H(:,n+1) = 0;                
            end
            
            % 6. Invert H to get the covariance matrix. H may be singular. 
            % Turn warning MATLAB:singularMatrix off, do the inversion and 
            % restore the warning state.
            warnState = warning('query','all');
            warning('off','MATLAB:singularMatrix');
            warning('off','MATLAB:illConditionedMatrix');
            warning('off','MATLAB:nearlySingularMatrix');
            cleanupObj = onCleanup(@() warning(warnState));            
            try
                C = covarianceOnNaturalScale(sglme,H);
            catch ME %#ok<NASGU>
                C = H \ eye(size(H));
            end
            
            % 7. C is square and of size length(beta) + length(eta) + 1.
                         lenbeta = length(beta);
                          leneta = length(eta);                        
                           ibeta = 1:lenbeta;
                    ietalogsigma = lenbeta + 1 : lenbeta + leneta + 1;            
                      covbetaHat = C(ibeta,ibeta); 
            covetaHatlogsigmaHat = C(ietalogsigma,ietalogsigma);            
            
        end % covBetaHatetaHatlogsigmaHatAsFunctionOfBetaThetaSigma.
        
        function [covbetaHat,covetaHatlogsigmaHat] = covBetaHatEtaHatLogSigmaHat(sglme)
%[covbetaHat,covetaHatlogsigmaHat] = covBetaHatEtaHatLogSigmaHat(sglme)
%   computes the joint covariance of [betaHat;etaHat;log(sigmaHat)] using
%   the joint Hessian of the negative Laplacian log likelihood with respect
%   to [beta;eta;log(sigma)] evaluated at [betaHat;etaHat;log(sigmaHat)].
%
%   If sglme.isSigmaFixed is true, then the joint covariance of
%   [betaHat;etaHat] is computed using the joint Hessian of the negative
%   Laplacian log likelihood with respect to [beta;eta] evaluated at
%   [betaHat;etaHat]. A row of all zeros and a column of all zeros is added
%   to this covariance to indicate that sigmaHat is estimated with 0
%   variance (i.e., it is fixed).            

            % 1. Get betaHat, thetaHat, sigmaHat.
             betaHat = sglme.betaHat;
            thetaHat = sglme.thetaHat;
            sigmaHat = sglme.sigmaHat;

            % 2. Get etaHat.
               Psi = sglme.Psi;
               Psi = setUnconstrainedParameters(Psi,thetaHat);
               Psi = setSigma(Psi,sigmaHat);
            etaHat = getNaturalParameters(Psi);
             
            % 3. If sigmaHat has been fixed then x = [betaHat;etaHat]
            % otherwise x = [betaHat;etaHat;log(sigmaHat)].                                                               
            if sglme.isSigmaFixed                
                x = [betaHat;etaHat];                
            else
                x = [betaHat;etaHat;log(sigmaHat)];
            end

            % 4. Make negative Laplacian log likelihood as a function of 
            % [beta;eta] if sigmaHat is fixed or [beta;eta;log(sigma)] if 
            % sigmaHat is not fixed.
            fun = makeNegativeLaplacianLogLikelihoodNaturalParameters(sglme);                       
            
            % 5. Get the Hessian of fun w.r.t [beta;eta] or w.r.t
            % [beta;eta;log(sigma)] evaluated at x.
            try
                H = sglme.getHessian(fun,x);
            catch ME %#ok<NASGU>
                H = NaN(length(x));
            end
            
            % 6. Invert H to get the covariance matrix. H may be singular. 
            % Turn warning MATLAB:singularMatrix off, do the inversion and 
            % restore the warning state.
            warnState = warning('query','all');
            warning('off','MATLAB:singularMatrix');
            warning('off','MATLAB:illConditionedMatrix');
            warning('off','MATLAB:nearlySingularMatrix');
            cleanupObj = onCleanup(@() warning(warnState));            
            try                
                C = sglme.covarianceOnNaturalScale(H);                
            catch ME %#ok<NASGU>
                C = H \ eye(size(H));
            end
            
            % 7. If C is the covariance of [betaHat;etaHat] then append a
            % column and row of all zeros to C.                        
            if sglme.isSigmaFixed    
                n = size(C,1);
                C(n+1,:) = 0;
                C(:,n+1) = 0;                
            end                        
            
            % 8. C is square and of size length(betaHat)+length(etaHat)+1.
                         lenbeta = length(betaHat);
                          leneta = length(etaHat);                        
                           ibeta = 1:lenbeta;
                    ietalogsigma = lenbeta + 1 : lenbeta + leneta + 1;            
                      covbetaHat = C(ibeta,ibeta); 
            covetaHatlogsigmaHat = C(ietalogsigma,ietalogsigma);            
            
        end % end of covBetaHatEtaHatLogSigmaHat.
                  
        function [pred,varpred] = getEstimateAndVariance(sglme,Xnew,Znew,betaBar,DeltabBar,theta,sigma,type)
%[pred,varpred] = getEstimateAndVariance(sglme,Xnew,Znew,betaBar,DeltabBar,theta,sigma)
%   takes specified values of (theta,sigma) and specified joint posterior 
%   mode approximation (betaBar,DeltabBar) corresponding to (theta,sigma)
%   and computes:
%
%   1. The posterior mode of each element of Xnew*beta + Znew*b given y,
%   theta and sigma in output pred.
%
%   2. The posterior variance of each element of Xnew*beta + Znew*b given
%   y, theta and sigma in output varpred.
%
%   Xnew is a M-by-p fixed effects design matrix. Znew is a M-by-q random 
%   effects design matrix. betaBar is p-by-1 and DeltabBar is q-by-1. theta
%   is the unconstrained parameter vector and sigma is square root of
%   dispersion parameter.
%
%   If M <= (p+q), interest may line in computing the posterior mode of the
%   "vector" Xnew*beta + Znew*b given y, theta and sigma and the associated
%   posterior covariance of the "vector" Xnew*beta + Znew*b given y, theta 
%   and sigma. If this is the intention, set the last input type to
%   'covariance'. 
%
%   Default value of type is 'variance' for computing the posterior
%   variance of each element of Xnew*beta + Znew*b given y, theta and
%   sigma.

            % 0. Default value of type.
            if nargin < 8
                type = 'variance';
            end

            % 1. Basic checks on inputs.
            M = size(Xnew,1);
            p = size(Xnew,2);
            q = size(Znew,2);
            assert( size(Znew,1) == M );            
            assert( p == sglme.p && q == sglme.q );
            assert( size(betaBar,1) == p && size(betaBar,2) == 1 );
            assert( size(DeltabBar,1) == q && size(DeltabBar,2) == 1 );

           % 2. Get X, N and Offset from sglme.
                X = sglme.X;                
                N = sglme.N;
            delta = sglme.Offset;
            
            % 3. Get U and Lambda matrices corresponding to theta.
            [U,Lambda] = getU(sglme,theta);
                        
            % 4. Get pred.
            pred = Xnew*betaBar + Znew*(Lambda*DeltabBar);
            pred = full(pred);
            
            % 5. Compute elementwise posterior variances if required.
            if nargout > 1
                warnState = warning('query','all');
                warning('off','MATLAB:nearlySingularMatrix');
                warning('off','MATLAB:illConditionedMatrix');
                warning('off','MATLAB:singularMatrix');
                cleanupObj = onCleanup(@() warning(warnState));
                
                % 5.1 Get the linear predictor using beta and Deltab.
                eta = X*betaBar + U*DeltabBar + delta;

                % 5.2 Compute the mean vector.
                ginv = sglme.Link.Inverse;
                  mu = ginv(eta);
                  mu = constrainMu(sglme,mu,sglme.Distribution);

                % 5.3. Get diagC and make it sparse diagonal matrix. For
                % FitMethod equal to 'mpl' or 'rempl', we replace the C
                % matrix by the W matrix to get the same results as from
                % the final fitted LME model from PL iterations.
                w = getEffectiveObservationWeights(sglme);
                switch lower(sglme.FitMethod)
                    case {'mpl','rempl'}
                        diagC = getDiagW(sglme,mu,w); 
                    otherwise
                        diagC = getDiagC(sglme,mu,w);
                end
                diagC = spdiags(diagC,0,N,N);   
                
                % 5.4 Get X'*C*X and X'*C*U
                XtCX = X'*(diagC*X);
                  CU = diagC*U;
                XtCU = X'*CU;

                % 5.5. Ensure that Cholesky factor of U'*C*U + eye(q)
                % exists. If not, the Laplace approximation is not well
                % defined and hence the joint covariance is also not well
                % defined.
                Iq = spdiags(ones(q,1),0,q,q);
                [R,status,S] = chol(U'*CU + Iq);
                if (status ~= 0)
                    error(message('stats:classreg:regr:lmeutils:StandardGeneralizedLinearMixedModel:BadLaplaceLogLikelihood'));
                end           
                
                % 5.6 For GLME models with non-canonical links, the C
                % matrix may have negative elements. This may make the
                % joint covariance matrix non positive definite. If this
                % happens, we error out.
                Q1 = (XtCU*S) / R;                
                [R1,status1] = chol(XtCX - Q1*Q1','lower');
                if (status1 ~= 0)
                    error(message('stats:classreg:regr:lmeutils:StandardGeneralizedLinearMixedModel:BadJointPosteriorCovariance'));
                end
                
                % 5.7 Compute Cb and Ca.
                if ( M <= q )
                    Cb = R' \ (S'*(Lambda'*Znew')); 
                else
                    Cb = (R' \ (S'*Lambda')) * Znew';                
                end
                Ca = R1 \ (Xnew' - Q1*Cb);
                
                % 5.8 Compute varpred and make it a column vector.
                if strcmpi(type,'variance')
                    varpred = sum(Cb.^2,1) + sum(Ca.^2,1);
                    varpred = (sigma^2) * varpred';      
                else
                    varpred = sigma^2*(Cb'*Cb + Ca'*Ca);
                end
                varpred = full(varpred);
            end
            
        end % end of getEstimateAndVariance.
        
    end
    
% Protected methods related to computing the covariance of 
% [thetaHat;log(sigmaHat)] and [etaHat,log(sigmaHat)]. 
    methods (Access=protected)
        
% NOTE: It is not required to compute this covariance for GLMEs. Hence the 
% following lines are commented out.        
%         function C = covThetaHatLogSigmaHat(slme)
%             
%         end % end of covThetaHatLogSigmaHat.
        
        function C = covEtaHatLogSigmaHat(sglme)
%C = covEtaHatLogSigmaHat(sglme) computes the covariance of [etaHat;log(sigmaHat)].            

            [~,C] = covBetaHatEtaHatLogSigmaHat(sglme);
            
        end % end of covEtaHatLogSigmaHat.
                
    end
    
% Public methods related to fitting.
    methods (Access=public)                
        
        function sglme = StandardGeneralizedLinearMixedModel(X,y,Z,Psi,FitMethod,dofit,dostats,varargin)
%

%StandardGeneralizedLinearMixedModel - Create a GLME model in standard form.
%   sglme = StandardGeneralizedLinearMixedModel(X,y,Z,Psi,FitMethod,dofit,dostats)
%
%           X = N-by-p fixed effects design matrix
%           y = N-by-1 response vector
%           Z = N-by-q random effects design matrix
%         Psi = an object of type CovarianceMatrix encapsulating the 
%               covariance matrix (sigma^2 * D)
%   FitMethod = a string that is one of 'MPL', 'REMPL', 'Laplace',
%               'ApproximateLaplace' or 'Quadrature'.
%
%   Returns a fitted StandardGeneralizedLinearMixedModel object sglme if
%   dofit is true. y is assumed to be conditionally independent given b. If
%   dofit is false then an unfitted object is returned which you must then
%   fit by calling the refit method. If dostats is true then the returned
%   object is ready for doing stats otherwise you must call the initstats
%   method on the object to make it ready for doing stats. If dofit is
%   false then an unfitted object is returned and dostats has no effect.
%
%   In most cases, you would set both dofit and dostats to true. If you
%   just want to access the estimated parameters such as thetaHat, betaHat,
%   bHat, sigmaHat, loglikHat, Psi etc. then set dofit to true and dostats
%   to false. The idea is that computing stats may not be needed for
%   certain operations such as likelihood ratio tests.
%
%   Before calling any of the stats methods, both isFitToData and
%   isReadyForStats flags must be true.
%
%   sglme = StandardGeneralizedLinearMixedModel(...,Name,Value,...) also
%   supplies optional name/value pairs:
%
%           Name                      Value
%           'Distribution'            Name of distribution for modelling 
%                                     the conditional distribution of y 
%                                     given b chosen from the following:
%                 'normal'              Normal distribution (default)
%                 'binomial'            Binomial distribution
%                 'poisson'             Poisson distribution
%                 'gamma'               Gamma distribution
%                 'inverse gaussian'    Inverse Gaussian distribution
%
%           'Link'                     The link function g to use in place
%                                      of the canonical link. Specify the
%                                      link as:
%                 'identity'    g(mu) = mu
%                 'log'         g(mu) = log(mu)
%                 'logit'       g(mu) = log(mu/(1-mu))
%                 'probit'      g(mu) = norminv(mu)
%                 'comploglog'  g(mu) = log(-log(1-mu))
%                 'loglog'      g(mu) = log(-log(mu))
%                 'reciprocal'  g(mu) = mu.^(-1)
%                 P (number)    g(mu) = mu.^P
%                 S (struct)    structure with four fields whose values
%                               are function handles and with these names:
%                                  S.Link               link function
%                                  S.Derivative         derivative
%                                  S.SecondDerivative   second derivative
%                                  S.Inverse            inverse of link
%                         The default is the canonical link that depends on
%                         the distribution:
%                 'identity'    normal distribution
%                 'logit'       binomial distribution
%                 'log'         Poisson distribution
%                 -1            gamma distribution
%                 -2            inverse gaussian distribution
%
%           'Offset'                   Vector of size N-by-1. This is used 
%                                      as an additional predictor with a
%                                      coefficient value fixed at 1.0.
%                                      Default is zeros(N,1).
%
%           'BinomialSize'             Vector of size N-by-1, specifying 
%                                      the size of the sample (number of
%                                      trials) used in computing y. This is
%                                      accepted when the 'Distribution'
%                                      parameter is 'binomial'. Default is
%                                      ones(N,1).
%
%           'Weights'                  Vector of size N-by-1. Default
%                                      is ones(N,1). For binomial and
%                                      Poisson distributions, 'Weights'
%                                      must be a vector of positive
%                                      integers.
%
%           'DispersionFlag'           Either true or false. Applies if 
%                                      'FitMethod' is 'MPL' or 'REMPL'. If
%                                      false, the dispersion parameter is
%                                      fixed at its theoretical value of
%                                      1.0 for binomial and Poisson
%                                      distributions. If true, the
%                                      dispersion parameter is estimated
%                                      from data even for binomial and
%                                      Poisson distributions. For all other
%                                      distributions, the dispersion
%                                      parameter is always estimated from
%                                      data. If 'FitMethod' is
%                                      'ApproximateLaplace', 'Laplace' or
%                                      'Quadrature', the dispersion
%                                      parameter is always fixed at 1.0 for
%                                      binomial and Poisson distributions
%                                      and estimated from data for all
%                                      other distributions. Default is
%                                      false.
%   
%           'Optimizer'               A string specifying the name of the 
%                                     optimizer. Supported values are
%                                     'quasinewton', 'fminunc' and
%                                     'fminsearch'. Default is
%                                     'quasinewton' if FitMethod is 'MPL'
%                                     or 'REMPL' and 'fminsearch'
%                                     otherwise.
%
%           'OptimizerOptions'        A structure containing the 
%                                     optimization options to be passed to
%                                     the optimizer. Created using
%                                     statset('fitglme') for 'quasinewton',
%                                     optimoptions('fminunc') for 'fminunc'
%                                     and optimset('fminsearch') for
%                                     'fminsearch'.
%
%           'InitializationMethod'    A string indicating how the initial
%                                     value of parameters should be chosen 
%                                     to initialize the optimizer. Either
%                                     'default' (default) or 'random'. 
%                                     Maximum likelihood based methods such
%                                     as 'Laplace', 'ApproximateLaplace'
%                                     and 'Quadrature' are always
%                                     initialized using 1 or more PL
%                                     iterations. 'InitializationMethod'
%                                     indicates the method used to
%                                     initialize the first PL iteration.
%
%           'CheckHessian'            Either true or false. If true, we
%                                     perform positive definiteness checks
%                                     on (a) Hessian of the objective
%                                     function at convergence (b) Covariance 
%                                     of [thetaHat;log(sigmaHat)] and (c)
%                                     Covariance of [etaHat;log(sigmaHat)].
%                                     If phiHat is also estimated then this
%                                     is also included in (b) and (c). 
%                                     Default is false.
%
%           'PLIterations'            A positive integer specifying the
%                                     maximum number of pseudo likelihood
%                                     (PL) iterations. Default is 100. PL
%                                     is used for fitting the model if
%                                     FitMethod is 'MPL' or 'REMPL'. For
%                                     other FitMethod values, PL iterations
%                                     are used to initialize parameters for
%                                     subsequent optimization.
%
%           'PLTolerance'             A real scalar to be used as a 
%                                     relative tolerance factor on the
%                                     linear predictor during PL 
%                                     iterations. Default is 1e-8.
%
%           'MuStart'                 A vector of size N-by-1 providing a
%                                     starting value for the conditional 
%                                     mean of y given b to initialize PL 
%                                     iterations. Legal values of elements 
%                                     in 'MuStart' are as follows:
%                 Distribution          Legal values
%                 'normal'              (-Inf,Inf)
%                 'binomial'            (0,1)
%                 'poisson'             (0,Inf)
%                 'gamma'               (0,Inf)
%                 'inverse gaussian'    (0,Inf)
%
%           'InitPLIterations'        Initial number of PL iterations used
%                                     to initialize parameters for maximum
%                                     likelihood (ML) based methods such as
%                                     'Laplace', 'ApproximateLaplace' and
%                                     'Quadrature'. Default is 10. Must be
%                                     greater than or equal to 1.
%
%           'EBMethod'                Method used to approximate the
%                                     empirical Bayes (EB) estimates of the
%                                     random effects. Choices are
%                                     'Default', 'LineSearchNewton',
%                                     'LineSearchModifiedNewton' and
%                                     'TrustRegion2D'. The 'Default'
%                                     'EBMethod' is similar to
%                                     'LineSearchNewton' but uses a
%                                     different convergence criterion and
%                                     does not display iterative progress.
%                                     'Default' and 'LineSearchNewton' may
%                                     fail for non-canonical link
%                                     functions. In these cases,
%                                     'TrustRegion2D' is the recommended
%                                     'EBMethod'.
%
%           'EBOptions'               A structure containing options for EB
%                                     optimization. The following options
%                                     are used:
%
%                 'TolFun'      Relative tolerance on the gradient norm.
%                               Default is 1e-6.
%                 'TolX'        Absolute tolerance on the step size.
%                               Default is 1e-8.
%                 'MaxIter'     Maximum number of iterations. Default is
%                               100.
%                 'Display'     'off', 'iter' or 'final'. Default is 'off'.
%
%                                     For 'EBMethod' equal to 'Default',
%                                     'TolFun' is the relative tolerance on
%                                     the linear predictor of the model and
%                                     the option 'Display' does not apply.
%
%           'CovarianceMethod'        Method used to compute the covariance
%                                     of estimated parameters. Choices are
%                                     'Conditional' (default) and
%                                     'JointHessian. The 'Conditional'
%                                     method computes a fast approximation
%                                     to the covariance of fixed effects
%                                     given the estimated covariance
%                                     parameters. The 'Conditional' method
%                                     does not compute the covariance of
%                                     covariance parameters. The
%                                     'JointHessian' method computes the
%                                     joint covariance of fixed effects and
%                                     covariance parameters via the
%                                     observed information matrix using the
%                                     Laplacian log likelihood.
%
%           'UseAMDPreordering'       Either true or false. If true, an
%                                     approximate minimum degree (AMD)
%                                     preordering is applied to the rows
%                                     and columns of the Jacobian used to
%                                     compute posterior modes of random
%                                     effects. Setting 'UseAMDPreordering'
%                                     to true can result in significant
%                                     speedup for larger models. Default is
%                                     false.
%
%           'NewtonStepMethod'        A string indicating the method to use 
%                                     for computing the Newton step during 
%                                     posterior mode estimation. Valid 
%                                     values are:
%
%                 'Cholesky'    Attempt Cholesky factorization of 
%                               (U'*C*U + I_q) first (possibly with AMD 
%                               preordering) and if not successful, use 
%                               backslash \. This is the default.
%                 'Backslash'   Always use backslash \.
%
%
%           'MuLowerBound'            Lower bound imposed on conditional
%                                     mean of y given b during iterations.
%                                     See the MuBound property for details.
%                                     
%           'MuUpperBound'            Upper bound imposed on conditional
%                                     mean of y given b during iterations.
%                                     See the MuBound property for details.
%
%           'UseSequentialFitting'    Either true or false. If false, all
%                                     ML methods are initialized using 1 or
%                                     more PL iterations. If true, the
%                                     initial values from PL iterations are
%                                     refined with 'ApproximateLaplace' for
%                                     'Laplace' based fitting and with
%                                     'ApproximateLaplace' followed by
%                                     'Laplace' for 'Quadrature' based
%                                     fitting. Default is false.
%
%           'ShowPLOptimizerDisplay'  Either true or false. If true, we
%                                     show optimizer display from PL
%                                     iterations. If false, we suppress
%                                     optimizer display from PL iterations.
%                                     Default is false.

            % (0) No arg constructor
            if ( nargin == 0 )
                return;
            end

            % (1) Ensure that dofit and dostats are scalar logicals.
            assert( isscalar(dofit) & islogical(dofit) );
            assert( isscalar(dostats) & islogical(dostats) );
            
            % (2) Validate inputs.
            [X,y,Z,Psi,FitMethod] = validateInputs(sglme,X,y,Z,Psi,FitMethod);
            
            % (3) Set N, p and q.
            [N,p] = size(X); %#ok<*PROP>
            q = size(Z,2);            
            sglme.N = N;
            sglme.p = p;
            sglme.q = q;
            
            % (4) Set X, y, Z, Psi and FitMethod.
            sglme.X = X;
            sglme.y = y;
            sglme.Z = Z;
            sglme.Psi = Psi;
            sglme.FitMethod = FitMethod;            
            
            % (5) Parse optional name/value pairs.            
                % 5(a) Define parameter defaults.
                dfltDistribution = 'normal';
                dfltLink = [];
                dfltOffset = zeros(N,1);
                dfltBinomialSize = ones(N,1);
                dfltWeights = ones(N,1);
                dfltDispersionFlag = false;
                switch lower(FitMethod)
                    case {'mpl','rempl'}
                        dfltOptimizer = 'quasinewton';
                    otherwise
                        dfltOptimizer = 'fminsearch';
                end
                dfltOptimizerOptions = struct([]);
                dfltInitializationMethod = 'default';
                dfltCheckHessian = false;    
                dfltPLIterations = 100;
                dfltPLTolerance = 1e-8;
                dfltMuStart = [];               
                dfltInitPLIterations = 10;
                dfltEBMethod = 'default';
                dfltEBOptions = statset('TolFun',1e-6,'TolX',1e-8,'MaxIter',100,'Display','off');                                
                dfltCovarianceMethod = 'conditional';                                
                dfltUseAMDPreordering = false;
                dfltNewtonStepMethod = 'cholesky';
                dfltMuLowerBound = eps;
                dfltMuUpperBound = Inf;
                dfltUseSequentialFitting = false;
                dfltShowPLOptimizerDisplay = false;
                
                % 5(b) Optional parameter names and their default values.
                paramNames = {  'Distribution',   'Link',   'Offset',   'BinomialSize',   'Weights',   'DispersionFlag',   'Optimizer',   'OptimizerOptions',   'InitializationMethod',   'CheckHessian',   'PLIterations',   'PLTolerance',   'MuStart',   'InitPLIterations',   'EBMethod',   'EBOptions',   'CovarianceMethod',   'UseAMDPreordering',   'NewtonStepMethod',   'MuLowerBound',   'MuUpperBound',   'UseSequentialFitting',   'ShowPLOptimizerDisplay'};
                paramDflts = {dfltDistribution, dfltLink, dfltOffset, dfltBinomialSize, dfltWeights, dfltDispersionFlag, dfltOptimizer, dfltOptimizerOptions, dfltInitializationMethod, dfltCheckHessian, dfltPLIterations, dfltPLTolerance, dfltMuStart, dfltInitPLIterations, dfltEBMethod, dfltEBOptions, dfltCovarianceMethod, dfltUseAMDPreordering, dfltNewtonStepMethod, dfltMuLowerBound, dfltMuUpperBound, dfltUseSequentialFitting, dfltShowPLOptimizerDisplay};
                
                % 5(c) Parse optional parameter name/value pairs.                
                [distribution,linkSpec,offset,...
                    binomialsize,weights,dispersionflag,...
                    optimizer,optimizeroptions,...
                    initializationmethod,checkhessian,...
                    pliterations,pltolerance,mustart,...                    
                    initpliterations,...
                    ebmethod, eboptions,...
                    covariancemethod,...
                    useamdpreordering,...
                    newtonstepmethod,...
                    mulowerbound,muupperbound,...
                    usesequentialfitting,...
                    showploptimizerdisplay] ...
                    = internal.stats.parseArgs(paramNames,paramDflts,varargin{:});                

            % (6) Validate optional inputs.                                               
                % 6(a) Validate distribution.
                distribution = internal.stats.getParamVal(distribution,sglme.AllowedDistributions,'Distribution');                                
                                
                % 6(b) Validate link specification. Ensure that selected
                % link makes sense for the chosen distribution.
                if isempty(linkSpec)
                    % No link specified, use a default based on the
                    % specified distribution.
                    linkSpec = defaultLink(sglme,distribution);
                end
                linkStruct = sglme.validateLink(linkSpec,sglme.FitMethod);
                
                checkDistributionLinkCombination(sglme,distribution,linkStruct);
                
                % 6(c) Validate offset.
                offset = validateOffset(sglme,offset,sglme.N);
                
                % 6(d) Validate binomialsize.
                if strcmpi(distribution,'binomial')
                    binomialsize = validateBinomialSize(sglme,binomialsize,sglme.N);
                end
                
                % 6(e) Validate prior weights.
                weights = validateWeights(sglme,weights,sglme.N,distribution);
                
                % 6(f) Validate dispersionflag.
                dispersionflag = internal.stats.parseOnOff(dispersionflag,'DispersionFlag');                                                                                                              
                
                % 6(g) Validate the values in y based on the specified
                % distribution, binomialsize and weights.                
                validateyRange(sglme,sglme.y,binomialsize,weights,distribution);                
                                
                % 6(h) validate optimization related options.
                [optimizer,optimizeroptions,initializationmethod] ...
                    = validateOptimizationOptions(sglme,optimizer,optimizeroptions,initializationmethod);                               
                                                
                % 6(i) Validate checkhessian.
                checkhessian = internal.stats.parseOnOff(checkhessian,'CheckHessian');
                
                % 6(j) Validate pliterations.
                pliterations = validatePLIterations(sglme,pliterations);
                
                % 6(k) Validate pltolerance.
                pltolerance = validatePLTolerance(sglme,pltolerance);
                
                % 6(l) Validate mustart.
                mustart = validateMuStart(sglme,mustart,distribution,N);
                
                % 6(m) Validate initial PL iterations.
                initpliterations = validateInitPLIterations(sglme,initpliterations);
                
                % 6(n) Validate eb iteration parameters.                
                [ebmethod,eboptions] = validateEBParameters(sglme,ebmethod,eboptions,dfltEBOptions);
                
                % 6(o) Validate covariancemethod.
                covariancemethod = internal.stats.getParamVal(covariancemethod,sglme.AllowedCovarianceMethods,'CovarianceMethod');                
                
                % 6(p) Validate useamdpreordering.
                useamdpreordering = internal.stats.parseOnOff(useamdpreordering,'UseAMDPreordering');
                
                % 6(q) Validate newtonstepmethod.
                newtonstepmethod = internal.stats.getParamVal(newtonstepmethod,sglme.AllowedNewtonStepMethods,'NewtonStepMethod');
                
                % 6(r) Validate mulowerbound and muupperbound.
                [mulowerbound,muupperbound] = validateMuBounds(sglme,mulowerbound,muupperbound);                
                
                % 6(s) Validate usesequentialfitting
                usesequentialfitting = internal.stats.parseOnOff(usesequentialfitting,'UseSequentialFitting');
                
                % 6(t) Validate showploptimizerdisplay.
                showploptimizerdisplay = internal.stats.parseOnOff(showploptimizerdisplay,'ShowPLOptimizerDisplay');
                
            % (7) Set selected options in the object.                        
            sglme.Distribution = distribution;
            sglme.Link = linkStruct;            
            sglme.Offset = offset;
            sglme.BinomialSize = binomialsize;                        
            sglme.PriorWeights = weights;              
            sglme.DispersionFixed = ...
                setDispersionFixed(sglme,dispersionflag,distribution,FitMethod);    
            if ( sglme.DispersionFixed == true )
                sglme.isSigmaFixed = true;
                sglme.sigmaFixed = 1.0;
            end
            
            sglme.Optimizer = optimizer;
            sglme.OptimizerOptions = optimizeroptions;
            sglme.InitializationMethod = initializationmethod;
            sglme.CheckHessian = checkhessian;
            sglme.PLIterations = pliterations;
            sglme.PLTolerance = pltolerance;
            sglme.MuStart = mustart;
            
            sglme.InitPLIterations = initpliterations;
            
            sglme.EBMethod = ebmethod;
            sglme.EBOptions = eboptions;            
            
            sglme.CovarianceMethod = covariancemethod;
            
            sglme.UseAMDPreordering = useamdpreordering;
            
            %       'NewtonStepMethod'          'NewtonStepMethodCode'
            %          'Cholesky'                        1
            %          'Backslash'                       2  
            sglme.NewtonStepMethod = newtonstepmethod;
            switch lower(newtonstepmethod)
                case 'cholesky'                    
                    sglme.NewtonStepMethodCode = 1;
                case 'backslash'
                    sglme.NewtonStepMethodCode = 2;                    
            end
            
            sglme.MuBound.TINY = mulowerbound;
            sglme.MuBound.BIG  = muupperbound;
            
            sglme.UseSequentialFitting = usesequentialfitting;
            
            sglme.ShowPLOptimizerDisplay = showploptimizerdisplay;
            
            % (8) Get variance function info into the object.
            sglme.VarianceFunction = varianceFunction(sglme,sglme.Distribution);
            
            % (9) Fit the standard LME model if asked.
            if (dofit == true)
                sglme = refit(sglme);
                % (10) Make the object ready for doing stats if required.
                if (dostats == true)
                    sglme = initstats(sglme);
                end
            end                   
            
        end % end of StandardGeneralizedLinearMixedModel.        
        
        function sglme = refit(sglme)
%refit - Fit and refit a Generalized Linear Mixed Effects (GLME) model in standard form.
%   sglme = refit(sglme) fits a standard GLME model previously initialized
%   with the StandardGeneralizedLinearMixedModel method. The public
%   properties of sglme such as y, X, etc. can be changed if desired and
%   the model can be refit using the refit method. When changing X, y, Z,
%   Psi or FitMethod, a call to refit and initstats is necessary before
%   calling any of the methods that compute stats on the fitted model.
%
%   When refitting a previously fitted model, parameters should be
%   initialized at the previously estimated values.

            % 1. Do either PL or ML.
            switch lower(sglme.FitMethod)                
                case {'mpl','rempl'}
                    % Fit using pseudo likelihood (PL).
                    numIter = sglme.PLIterations;
                    kappa = sglme.PLTolerance;
                    [sglme,cause] = fitUsingPL(sglme,numIter,kappa);                    
                    if ( cause == 1 )
                        warning(message('stats:classreg:regr:lmeutils:StandardGeneralizedLinearMixedModel:Message_PLUnableToConverge',numIter));                           
                    end                    
                    
                case {'approximatelaplace','laplace','quadrature'}
                    % Fit using maximum likelihood (ML).
                    [sglme,cause] = fitUsingML(sglme);
                    if ( cause ~= 0 && cause ~= 1 )
                        warning(message('stats:classreg:regr:lmeutils:StandardGeneralizedLinearMixedModel:Message_OptimizerUnableToConverge',sglme.Optimizer));
                    end
                    
                otherwise
                    error(message('stats:classreg:regr:lmeutils:StandardGeneralizedLinearMixedModel:BadFitMethod'));                
            end                                    
            
            % 2. Mark the object as Fitted.
            sglme.isFitToData = true;
            
        end % end of refit.        
                
        function sglme = initstats(sglme)
%initstats - Make a fitted StandardGeneralizedLinearMixedModel object ready for stats.
%   sglme = initstats(sglme) takes a fitted GLME object sglme and makes it
%   ready for stats. sglme must be a fitted object returned by either
%   StandardGeneralizedLinearMixedModel constructor or the refit method.
%   This method is responsible for the following tasks:
%
%   1. Fill in covbetaHat.
%   2. Fill in covthetaHatlogsigmaHat.
%   3. Fill in covetaHatlogsigmaHat.
%   4. Do positive definiteness checks on covthetaHatlogsigmaHat and 
%      covetaHatlogsigmaHat if CheckHessian is true.
%   5. This method *does not* fill in covbetaHatbHat the covariance needed
%      for joint inference on beta/b. This is because covbetaHatbHat can be
%      a large matrix and so it makes sense to access it via a method on an
%      as needed basis.
%   6. Set isReadyForStats to true.

            % 1. Ensure that sglme has been fit to data.
            if (sglme.isFitToData == false)
                error(message('stats:classreg:regr:lmeutils:StandardGeneralizedLinearMixedModel:MustRefitFirst'));
            end

            % 2. initstatsPL simply calls the initstats method on the final
            % LME object from PL that is stored in sglme. initstatsML
            % follows a ML specific recipe.
            switch lower(sglme.FitMethod)
                case {'mpl','rempl'}
                    sglme = initstatsPL(sglme);
                case {'approximatelaplace','laplace','quadrature'}
                    sglme = initstatsML(sglme);
                otherwise
                    error(message('stats:classreg:regr:lmeutils:StandardGeneralizedLinearMixedModel:BadFitMethod'));
            end
            
            % 3. Mark the object as ready for stats.
            sglme.isReadyForStats = true;             

        end % end of initstats.                               
        
    end       
  
% Public methods to get Satterthwaite DF for T/F tests on beta/b vector.
% The following four methods simply return an error message since
% Satterthwaite DF computation for GLMEs may not be sensible.
    methods (Access=public)    
        
        function df = dfBetaTTest(slme,c) %#ok<STOUT,INUSD>
            
            error(message('stats:classreg:regr:lmeutils:StandardGeneralizedLinearMixedModel:BadDFMethod'));
            
        end
        
        function df = dfBetaFTest(slme,L) %#ok<STOUT,INUSD>
            
            error(message('stats:classreg:regr:lmeutils:StandardGeneralizedLinearMixedModel:BadDFMethod'));
            
        end
        
        function df = dfBetaBTTest(slme,c) %#ok<STOUT,INUSD>
            
            error(message('stats:classreg:regr:lmeutils:StandardGeneralizedLinearMixedModel:BadDFMethod'));
            
        end
        
        function df = dfBetaBFTest(slme,L) %#ok<STOUT,INUSD>
            
            error(message('stats:classreg:regr:lmeutils:StandardGeneralizedLinearMixedModel:BadDFMethod'));
            
        end        
        
    end
        
% Public method for computing model criterion.
    methods (Access=public)
       
        function crittable= modelCriterion(sglme)
%modelCriterion - Compute table containing model criterion info.
%   crittable= modelCriterion(sglme) takes a StandardGeneralizedLinearMixedModel object
%   sglme and computes a table of model criterion such as AIC, BIC, logLik 
%   and Deviance.

    % Create a structure stats for classreg.regr.modelutils.modelcriterion
    % such that:
    % L = stats.LogLikelihood;
    % k = stats.NumCoefficients;
    % n = stats.NumObservations;      
    
            % (1) Get N and p.
            N = sglme.N;
            p = sglme.p;
            
            % (2) Add the fixed effects + residual variance to
            % NumCoefficients.
            if ( sglme.isSigmaFixed == true )
                stats.NumCoefficients = ...
                sglme.Psi.NumParametersExcludingSigma + p;
            else
                stats.NumCoefficients = ...
                    sglme.Psi.NumParametersExcludingSigma + (p+1);
            end
                        
            % (3) Get effective number of observations based on FitMethod.
            switch lower(sglme.FitMethod)
                case {'mpl','approximatelaplace','laplace','quadrature'}
                    stats.NumObservations = N;
                case {'rempl'}
                    stats.NumObservations = (N-p);
                    
                otherwise
                    % <entry key="BadFitMethod">''FitMethod'' must be either ''ML'' or ''REML''.</entry>
                    error(message('stats:classreg:regr:lmeutils:StandardGeneralizedLinearMixedModel:BadFitMethod'));                    
            end
                                    
            % (4) Set maximized log likelihood.
            switch lower(sglme.FitMethod)
                case {'approximatelaplace','laplace','quadrature'}
                    stats.LogLikelihood = sglme.loglikHat;
                case {'mpl','rempl'}
                    stats.LogLikelihood = sglme.loglikHatPseudoData;                    
                otherwise
                    % <entry key="BadFitMethod">''FitMethod'' must be either ''ML'' or ''REML''.</entry>
                    error(message('stats:classreg:regr:lmeutils:StandardGeneralizedLinearMixedModel:BadFitMethod'));                    
            end
                                    
            % (5) Call modelcriterion utility function.
            crit = classreg.regr.modelutils.modelcriterion(stats,'all',true);
            
            % (6) Get deviance and create output table.
            Deviance = -2*stats.LogLikelihood;
            crittable = table(crit.AIC, crit.BIC, stats.LogLikelihood, Deviance,...
                'VariableNames',{'AIC' 'BIC' 'logLik' 'Deviance'});
            
        end % end of modelCriterion.      
        
    end    
        
% Public method to get estimated posterior covariance of random effects vector.
    methods (Access=public)
        
        function H = postCovb(sglme)
%postCovb - Posterior covariance of b. 
%   H = postCovb(sglme) takes a StandardGeneralizedLinearMixedModel object
%   sglme and computes the estimated posterior covariance matrix of b
%   conditional on y and estimated values of beta, theta, sigma^2.

            % The posterior covariance is given by:
            % H = sigma^2 * inv(Z'*C*Z + Delta'*Delta) 
            % with Delta'*Delta = inv(D).
            
            % 1. Get X and Offset from sglme.
                X = sglme.X;               
            delta = sglme.Offset;
            
            % 2. Get estimated model parameters.
             thetaHat = sglme.thetaHat;
              betaHat = sglme.betaHat;
            DeltabHat = sglme.DeltabHat;
             sigmaHat = sglme.sigmaHat;
            
            % 3. Get U and Lambda matrices at thetaHat.
            [U,Lambda] = getU(sglme,thetaHat);
                        
            % 4. Get the linear predictor using betaHat and DeltabHat.
            eta = X*betaHat + U*DeltabHat + delta;
            
            % 5. Compute the mean vector.
            ginv = sglme.Link.Inverse;
              mu = ginv(eta);
              mu = constrainMu(sglme,mu,sglme.Distribution);
              
            % 6. Get diagC. For FitMethod equal to 'mpl' or 'rempl', we
            % replace the C matrix by the W matrix to get the same results
            % as from the final fitted LME model from PL iterations.
            w = getEffectiveObservationWeights(sglme);
            switch lower(sglme.FitMethod)
                case {'mpl','rempl'}
                    diagC = getDiagW(sglme,mu,w); 
                otherwise
                    diagC = getDiagC(sglme,mu,w);
            end            
                        
            % 7. Get R and S such that S*R'*R*S' = U'*C*U + eye(q) where R
            % and S are q-by-q. S is a permutation matrix: S*S' = eye(q).                      
            M = getUtCUPlusIdentity(sglme,U,diagC);
            [R,status,S] = chol(M);
            
            % 8. Make sure status is zero.
            if (status ~= 0)
                error(message('stats:classreg:regr:lmeutils:StandardGeneralizedLinearMixedModel:BadLaplaceLogLikelihood'));
            end            
            
            % 9. H can also be written as:
            %    H = sigma^2 * (Lambda*S*inv(R) * inv(R)'*S'*Lambda').
            T = (Lambda*S) / R;
            H = (sigmaHat^2) * (T*T');
            
        end % end of postCovb.
    
    end        
  
% Public methods to get fitted values.
    methods (Access=public)   
        
        function mufit = fitted(sglme,wantConditional)
%fitted - Returns the fitted response from a StandardGeneralizedLinearMixedModel.
%   mufit = fitted(sglme,wantConditional) returns a N-by-1 vector mufit
%   that represents the fitted mean of the response in the model. If
%   wantConditional is true, mufit represents the fitted conditional mean
%   of the response. On the other hand, if wantConditional is false, mufit
%   represents the fitted conditional mean of the response with the
%   empirical Bayes predictor (EBP) of the random effects set to zero.

            % 1. Ensure that wantConditional is sensible.
            assert(islogical(wantConditional) ...
                & isscalar(wantConditional));

            % 2. Get X, betaHat, delta and inverse link.
                  X = sglme.X;
            betaHat = sglme.betaHat;
              delta = sglme.Offset;
               ginv = sglme.Link.Inverse;
            
            % 3. Form the appropriate linear predictor etaHat.
            if (wantConditional == true)
                % 3.1 Get Z and bHat.
                   Z = sglme.Z;
                bHat = sglme.bHat;
                % 3.2 etaHat using both betaHat and bHat.
                etaHat = X*betaHat + Z*bHat + delta;                
            else
                % 3.3 etaHat using only betaHat.
                etaHat = X*betaHat + delta;
            end
            
            % 4. Apply ginv to convert etaHat to the scale of response.
            mufit = ginv(etaHat);            
            mufit = constrainMu(sglme,mufit,sglme.Distribution);
            
        end % end of fitted.               
                               
    end
            
% Public method to make predictions on new data.
    methods (Access=public)
       
        function [ypred,CI,DF] = predict(sglme,Xnew,Znew,alpha,dfmethod,wantConditional,wantPointwise,delta,hasIntercept)
%predict - Predictions from a fitted StandardGeneralizedLinearMixedModel.
%   [ypred,CI,DF] = predict(sglme,Xnew,Znew,alpha,dfmethod,wantConditional,wantPointwise,delta,hasIntercept)
%   takes a fitted StandardGeneralizedLinearMixedModel sglme, a fixed
%   effects design matrix Xnew, a random effects design matrix Znew, a
%   confidence level alpha, a degrees of freedom method dfmethod, flags for
%   conditional and pointwise predictions, an offset vector delta and
%   outputs predictions and (1-alpha) confidence intervals (CIs) for the
%   true predictions from sglme at Xnew and Znew. If wantConditional is
%   true then the predictions include contributions from both fixed effects
%   and EBPs of random effects otherwise the predictions include
%   contributions from only the fixed effects. If wantPointwise is true
%   then pointwise CIs are computed otherwise simultaneous CIs are
%   computed. Xnew must be M-by-p and Znew must be M-by-q where p = sglme.p
%   and q = sglme.q. The offset vector delta is M-by-1. The output ypred is
%   a M-by-1 vector containing the predictions corresponding to the rows of
%   Xnew and Znew and output CI is a M-by-2 matrix. Each row of CI is a
%   (1-alpha) confidence interval for the true prediction such that CI(:,1)
%   contains the lower confidence limit and CI(:,2) contains the upper
%   confidence limit. DF is a M-by-1 vector containing the DF values used
%   in computing the CIs if wantPointwise is true. If wantPointwise is
%   false then DF is a scalar containing the DF value used for computing
%   the Scheffe simultaneous CIs. hasIntercept should be set to true if the
%   fixed effects part of the model has an intercept and false otherwise.
            
            % 1. Enforce that we never return 'observation' predictions for
            % a GLME and get predictions and CIs on the GLME linear
            % predictor.
            wantCurve = true;
            args = {Xnew,Znew,alpha,dfmethod,wantConditional,wantPointwise,wantCurve,hasIntercept};
            switch nargout
                case {0,1}
                    ypred         = predict@classreg.regr.lmeutils.StandardLinearLikeMixedModel(sglme,args{:});
                case 2
                    [ypred,CI]    = predict@classreg.regr.lmeutils.StandardLinearLikeMixedModel(sglme,args{:});
                case 3
                    [ypred,CI,DF] = predict@classreg.regr.lmeutils.StandardLinearLikeMixedModel(sglme,args{:});
            end
            %[ypred,CI,DF] = predict@classreg.regr.lmeutils.StandardLinearLikeMixedModel(sglme,...
            %    Xnew,Znew,alpha,dfmethod,wantConditional,wantPointwise,wantCurve,hasIntercept);
            
            % 2. Validate the offset vector delta.
            assert( isnumeric(delta) & isreal(delta) & iscolumn(delta) & size(delta,1) == size(ypred,1) );
            
            % 3. Add the offset to ypred and CI and transform the CIs to
            % the scale of the modelled response.
                    ginv = sglme.Link.Inverse;        
            distribution = sglme.Distribution;
                   ypred = ginv(ypred + delta); 
                   ypred = constrainMu(sglme,ypred,distribution);
            if nargout > 1
                     CI = ginv(bsxfun(@plus,CI,delta));     
                CI(:,1) = constrainMu(sglme,CI(:,1),distribution);
                CI(:,2) = constrainMu(sglme,CI(:,2),distribution);                
                     CI = [min(CI,[],2),max(CI,[],2)];
            end
            
        end % end of predict.
                                
    end
    
% Public method for generating random data from the fitted model.
    methods (Access=public)
                
        function ysim = random(sglme,S,Xsim,Zsim,delta,wp,ntrials,numreps)
%random - Generate random data from a fitted StandardGeneralizedLinearMixedModel.
%   ysim = random(sglme,S,Xsim,Zsim,delta,wp,ntrials) generates random data
%   from a fitted StandardGeneralizedLinearMixedModel object sglme using
%   the RandStream object S at the specified fixed effects design matrix
%   Xsim, random effects design matrix Zsim, offset vector delta, prior
%   weights vector wp and for 'binomial' distribution the number of trials
%   per observation ntrials. Xsim must be M-by-p and Zsim must be M-by-q
%   where p = sglme.p and q = sglme.q. Vectors delta, wp and ntrials are of
%   sizes M-by-1 each. The output ysim is a M-by-1 vector. Set S = [] to 
%   use the default RandStream object. For non 'binomial' distributions,
%   ntrials may be set []. Arguments delta, wp and ntrials are optional
%   with the default values of zeros(M,1), ones(M,1) and ones(M,1)
%   respectively.
%
%   ysim = random(sglme,S,Xsim,Zsim,delta,wp,ntrials,numreps) also takes a 
%   NumBlocks-by-1 vector numreps where NumBlocks = sglme.Psi.NumBlocks. In 
%   this case, we assume that Zsim is of size:
%
%   size(Zsim,2) = sum_{i=1 to NumBlocks} SizeVec(i) * numreps(i)
%   
%   where SizeVec = sglme.Psi.SizeVec. The random effects vector bsim
%   corresponding to Zsim is generated by the call:
%
%   bsim = randomb(sglme,S,numreps)
%
%   If numreps is not supplied, it is taken to be sglme.Psi.NumReps and in
%   that case Zsim must be of size M-by-q where q = sglme.q.

            % 1. Can be called with 4,5,6,7 or 8 input args.
            narginchk(4,8);

            % 2. Ensure that Xsim is sensible.
            assert( isnumeric(Xsim) & isreal(Xsim) & ismatrix(Xsim) );
            assert( size(Xsim,2) == sglme.p );
            
            % 3. Ensure that Zsim is sensible.
            assert( isnumeric(Zsim) & isreal(Zsim) & ismatrix(Zsim) );            
            
            % 4. Ensure that Xsim and Zsim have the same number of rows.
            assert( size(Xsim,1) == size(Zsim,1) );
            
            % 5. Set default values for delta, wp and ntrials.
            N = size(Xsim,1); 
            if nargin < 5
                delta = zeros(N,1);
            end
            if nargin < 6
                wp = ones(N,1);
            end
            if nargin < 7
                ntrials = ones(N,1);
            end
            if nargin < 8
                numreps = sglme.Psi.NumReps;
            end
            
            % 6. Ensure delta is sensible.
            assert( isnumeric(delta) & isreal(delta) & iscolumn(delta) );
            
            % 7. Simulate vector bsim from sglme.Psi. The combined random
            % effects vector bsim will be of size sum_{i=1 to NumBlocks}
            % SizeVec(i) * numreps(i) where NumBlocks = sglme.Psi.NumBlocks
            % and SizeVec = sglme.Psi.SizeVec. If numreps is [] then bsim 
            % is a 0-by-1 vector.
            if isempty(numreps)
                bsim = zeros(0,1);
            else
                bsim = randomb(sglme,S,numreps);
            end
            assert( size(Zsim,2) == length(bsim) );
            
            % 8. Get betaHat and sigmaHat from sglme.                                    
             betaHat = sglme.betaHat;
            sigmaHat = sglme.sigmaHat;
                        
            % 9. Get the linear predictor eta and transform it to mu.
            % Ensure that mu satisfies distribution specific constraints.
            if isempty(bsim)
                eta = Xsim*betaHat + delta;
            else
                eta = Xsim*betaHat + Zsim*bsim + delta;
            end
            mu = sglme.Link.Inverse(eta);
            mu = constrainMu(sglme,mu,sglme.Distribution);
            
            % 10. Generate data from conditional distribution of y | bsim.
            % wp and ntrials are validated inside conditionalRandom.
            if isempty(S)
                 ysim = conditionalRandom(sglme,mu,sigmaHat,wp,ntrials);
            else
                prevS = RandStream.setGlobalStream(S);
           cleanupObj = onCleanup(@() RandStream.setGlobalStream(prevS));
                 ysim = conditionalRandom(sglme,mu,sigmaHat,wp,ntrials);                
            end                            

        end % end of random.
        
    end    

% Protected method to generate random numbers from the conditional 
% distribution of y | b using the mu, sigma parameterization.   
    methods (Access=protected)
       
        function yrnd = conditionalRandom(sglme,mu,sigma,wp,ntrials)
%yrnd = conditionalRandom(sglme,mu,sigma,wp,ntrials) generates random 
%   values from the distribution of y given b parameterized by M-by-1
%   conditional mean vector mu, square root of dispersion parameter sigma,
%   M-by-1 prior weight vector wp and for the 'binomial' distribution a
%   M-by-1 vector ntrials containing the number of trials per new
%   observation. For non 'binomial' distributions, ntrials can be set [].
            
            % 1. Get the current distribution.
            distribution = sglme.Distribution;

            % 2. Validate mu, sigma, wp.
            M = size(mu,1);
            mu =  validateMuStart(sglme,mu,distribution,M);
            sigma = validateSigma(sglme,sigma);
            wp =  validateWeights(sglme,wp,M,distribution);            
                                    
            % 3. Call random as required.
            switch lower(distribution)              
                case 'binomial'  
                    ntrials = validateBinomialSize(sglme,ntrials,M);                    
                    counts = ntrials.*wp;                                 
                    sumurnd = random('binomial', counts, mu, M, 1);
                    yrnd = sumurnd./counts;  
                case 'poisson'                    
                    sumurnd = random('poisson', wp.*mu, M, 1);
                    yrnd = sumurnd./wp;                    
                case 'gamma'
                    a = (1/(sigma^2));
                    awp = a.*wp;
                    b = mu./awp;
                    yrnd = random('gamma', awp, b, M, 1);
                case {'inverse gaussian','inversegaussian'}
                    lambda = (1/(sigma^2));
                    yrnd = random('inversegaussian', mu, lambda.*wp, M, 1);    
                case {'normal','gaussian'}
                    yrnd = random('normal', mu, sigma./sqrt(wp), M, 1);
            end

        end % end of conditionalRandom.
        
    end

% Protected methods to get raw/Pearson/Standardized residuals.     
    methods (Access=protected)        
        
        function rawr = getRawResiduals(sglme,wantConditional)
%getRawResiduals - Get Raw residuals.
%   rawr = getRawResiduals(sglme,wantConditional) takes an object sglme of
%   type StandardGeneralizedLinearMixedModel and returns the raw
%   conditional residuals if wantConditional is true and the raw marginal
%   residuals if the flag wantConditional is false.
           
            % 1. Ensure that wantConditional is sensible.
            assert(isscalar(wantConditional) & islogical(wantConditional));
            
            % 2. Form the raw residuals by subtracting the appropriate type
            % of fitted response from the response vector.
            mufit = fitted(sglme,wantConditional);
             rawr = sglme.y - mufit;
           
        end % end of getRawResiduals.
        
        function pearsonr = getPearsonResiduals(sglme,wantConditional)
%getPearsonResiduals - Get Pearson residuals.           
%   pearsonr = getPearsonResiduals(sglme,wantConditional) takes an object
%   sglme of type StandardGeneralizedLinearMixedModel and returns the
%   Pearson conditional residuals if wantConditional is true and the
%   Pearson marginal residuals if wantConditional is false.

            % 1. Ensure that wantConditional is sensible.
            assert(isscalar(wantConditional) ...
                & islogical(wantConditional));

            % 2. Get fitted values of the right type.
            mufit = fitted(sglme,wantConditional);
            
            % 3. Get the raw residuals based on mufit.
            rawr = sglme.y - mufit;
            
            % 4. Get the normalizing factor for rawr.
            sigmaHat = sglme.sigmaHat;
                 w = getEffectiveObservationWeights(sglme);
                 v = sglme.VarianceFunction.VarianceFunction;            
            stdvec = sigmaHat * sqrt(v(mufit) ./ w);
            
            % 5. Get the Pearson residuals.
            pearsonr = rawr ./ stdvec;                                   
            
        end % end of getPearsonResiduals.                   
        
        function studr = getStudentizedResiduals(sglme,wantConditional) %#ok<STOUT,INUSD>
                     
             % It is not clear how to define "Studentized" or
             % "Standardized" residuals for GLME models. Hence, this method
             % simply returns an error.
             error(message('stats:classreg:regr:lmeutils:StandardGeneralizedLinearMixedModel:BadResidualType'));
             
        end
        
    end
        
% Static utility methods.
    methods (Static)
        
        function optsFsolve = convertOptionsToFSolveOptions(opts,dfltopts)
%optsFsolve = convertOptionsToFSolveOptions(opts,dfltopts) takes a 
%   "partially filled" struct opts and a "fully filled" statset struct
%   dfltopts and returns a optiomoptions('fsolve') object containing filled
%   in values for 'TolFun', 'TolX', 'MaxIter' and 'Display'.

            % 1. Ensure that both inputs are structs.
            assert( isstruct(opts) & isstruct(dfltopts) );

            % 2. Combine opts and dfltopts to create a fully filled opts.
            opts = statset(dfltopts,opts);
                
            % 3. Create the default optimoptions object for fsolve.
            optsFsolve = optimoptions('fsolve'); 
            
            % 4. Copy the fields, TolFun, TolX, MaxIter and Display from
            % opts to optsFsolve.            
            optsFsolve.TolFun  = opts.TolFun;
            optsFsolve.TolX    = opts.TolX;
            optsFsolve.MaxIter = opts.MaxIter;
            optsFsolve.Display = opts.Display;  
            
        end % end of convertOptionsToFSolveOptions.
        
        function pn = computeNewtonStepForB(J,r,useamd,s,methodcode)
%pn = computeNewtonStep(J,r,useamd,s,methodcode) takes a q-by-q potentially
%   symmetric, positive definite matrix J, a q-by-1 vector r, a logical
%   flag useamd, a q-by-1 integer vector s and an integer methodcode and
%   solves J*pn = -r for pn. methodcode = 1 corresponds to 'Cholesky' and
%   methodcode = 2 corresponds to 'Backslash'. If useamd is false then s is
%   not used (and hence can be set to []).


            switch methodcode
                case 1
                    % 1. Try Cholesky and if that doesn't work use \.
                    if (useamd == true)
                        [L,status] = chol(J(s,s),'lower');
                    else
                        [L,status,s] = chol(J,'lower','vector');
                    end
                    if (status == 0)
                        % 2. J is positive definite.
                        pn = r;
                        pn(s) = -(L'\(L\r(s)));
                    else
                        % 3. Must use \. Use (-J)\r so that \ does not
                        % try Cholesky.
                        pn = (-J)\r;
                    end
                case 2
                    % 4. Always use \. Use (J\r) so that \ tries
                    % Cholesky.
                    pn = -(J\r);
                otherwise
            end

        end % end of computeNewtonStepForB.
        
        function pn = computeNewtonStepForBetaB(J,r,UtCUIq,XtCU,XtCX,rb,rbeta,p,q,Iq,useamd,s,methodcode)
%pn = computeNewtonStepForBetaB(J,r,UtCUIq,XtCU,XtCX,rb,rbeta,p,q,Iq,useamd,s,methodcode)
%   solves the system of equations for computing the Newton step in the
%   approximate Laplace method. J is the (p+q)-by-(p+q) Jacobian
%   corresponding to the (p+q)-by-1 system of non linear equations in r.
%   UtCUIq is U'*C*U + I_q, XtCU is X'*C*U, XtCX is X'*C*X, rb and rbeta
%   form the (p+q)-by-1 system of non-linear equations in the approximate
%   Laplace method. p is the number of fixed effects and q is the number of
%   random effects. Iq is the q-by-q identity matrix. useamd is either true
%   or false. methodcode is an integer such that methodcode = 1 corresponds
%   to 'Cholesky' and methodcode = 2 corresponds to 'Backslash'. s is a
%   q-by-1 integer vector. If useamd is false then s is not used (hence can
%   be set []).

            switch methodcode
                case 1   
                    % 1. Try to use Cholesky and if that doesn't work use
                    % \.
                    if (useamd == true)
                        [L,status] = chol(UtCUIq(s,s),'lower');
                    else
                        [L,status,s] = chol(UtCUIq,'lower','vector');
                    end                    
                    if (status == 0)    
                        % 2. UtCUIq is positive definite.
                                  S = Iq(:,s);
                                 Q1 = (XtCU*S)/L';
                                 cb = L\(S'*rb);
                                 pn = r;
                            pn(1:p) = (XtCX - Q1*Q1')\(rbeta - Q1*cb);
                        pn(p+1:p+q) = S*(L'\(cb - Q1'*pn(1:p)));
                                 pn = -pn;
                    else
                        % 3. Must use \.
                        pn = (-J)\r;
                    end                    
                case 2                  
                    % 4. Always use \.
                    pn = -(J\r);
                otherwise   
            end

        end % end of computeNewtonStepForBetaB.
        
    end
    
end
