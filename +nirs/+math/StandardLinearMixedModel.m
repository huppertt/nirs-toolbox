classdef StandardLinearMixedModel < classreg.regr.lmeutils.StandardLinearLikeMixedModel
%  This feature is intended for internal use only and is subject to change
%  at any time without warning. 

%StandardLinearMixedModel - An internal utility to fit Linear Mixed Effects
%   (LME) models in standard form. This feature is intended for internal 
%   use only and is subject to change at any time without warning. 
%%
% $y = X \beta + Z b + \varepsilon$
% 
% * $y$ is a $N \times 1$ vector.
% * $X$ is a $N \times p$ fixed effects design matrix
% * $\beta$ is a $p \times 1$ fixed effect vector
% * $Z$ is a $N \times q$ random effects design matrix
% * $b$ is a $q \times 1$ random effect vector
% * $\varepsilon \sim N(0,\sigma^2 I_N)$
% * $b \sim N(0, \sigma^2 D)$
% * $\Psi = \sigma^2 D$
%
%   The purpose of this class is to:
%   (1) Fit LME models in standard form by ML or REML.
%   (2) Refit the LME model after changing some inputs.
%   (3) Compute the estimate betaHat of beta.
%   (4) Compute the estimate thetaHat of theta.
%   (5) Compute the estimate sigmaHat of sigma.
%   (6) Compute the estimated BLUP bHat of b.
%   (7) Estimate covariance of betaHat.
%   (8) Perform hypothesis tests on beta.
%   (9) CIs for contrasts of beta.
%  (10) Compute covariance of [betaHat - beta;bHat - b].
%  (11) Peform hypothesis tests on [beta;b].
%  (12) CIs for contrasts of [beta;b].
%  (13) Estimate covariance of [thetaHat;log(sigmaHat)].
%  (14) Estimate covariance of [etaHat;log(sigmaHat)] where etaHat is the
%       Natural parameter vector for [thetaHat;log(sigmaHat)].
%  (15) Compute confidence intervals on the canonical parameter vector
%       [heta;log(sigma)].
%  (16) Compute conditional/marginal fitted values.
%  (17) Compute conditional/marginal residuals of various types.
%  (18) Generate random data from fitted model.
%  (19) Make predictions on new data.
%  (20) Compute the posterior covariance of random effects vector.
%  (21) Compute various goodness of fit criterion.
%
%   StandardLinearMixedModel methods:
%      StandardLinearMixedModel - Create/Fit a LME model in standard form.
%      refit - Fit/refit a previously initialized LME model.
%      initstats - Make a fitted LME model ready for stats.
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
%   StandardLinearMixedModel properties:
%      y - Response vector.
%      X - Fixed effects design matrix.
%      Z - Random effects design matrix.
%      Psi - An object of type CovarianceMatrix encapsulating $\sigma^2 D$.
%      FitMethod - Method used to fit the LME.  
%      N - Number of rows in y, X and Z.
%      p - Number of columns in X.       
%      q - Number of columns in Z.
%      rankX - Rank of the fixed effects design matrix X.
%      Optimizer - Selected optimizer.
%      OptimizerOptions - Options for the selected optimizer.
%      CheckHessian - Flag for Hessian checks at convergence.
%      InitializationMethod - Selected method for initializing parameters.
%      betaHat - Estimated fixed effects vector.        
%      bHat - Estimated Best Linear Unbiased Predictor of random effects.
%      DeltabHat - Estimated normalized posterior mode of random effects.
%      sigmaHat - Estimated residual standard deviation.       
%      thetaHat - Estimated unconstrained parameter for matrix D.
%      loglikHat - Maximized log likelihood or restricted log likelihood.   
%      covbetaHat - Estimated covariance of betaHat.
%      covthetaHatlogsigmaHat - Estimated covariance on unconstrained scale.       
%      covetaHatlogsigmaHat - Estimated covariance on Natural scale.        
%      covbetaHatbHat - Estimated covariance of [betaHat-beta;bHat-b].
%      isSigmaFixed - True if residual standard deviation sigma is fixed.
%      sigmaFixed - The specified fixed value of sigma.
%      isFitToData - Flag to mark a fitted StandardLinearMixedModel object.
%      isReadyForStats - Flag to mark the object ready for doing stats.

%   Copyright 2012-2013 The MathWorks, Inc.    

% Public properties supplied by the user.
    properties (GetAccess=public, SetAccess=public)
%y - N by 1 response vector used to fit the LME.        
        y

%X - N by p fixed effects design matrix used to fit the LME.        
        X

%Z - N by q random effects design matrix used to fit the LME.        
        Z

%Psi - An object of type CovarianceMatrix representing the matrix Psi =
%   sigma^2 * D(theta). If not using random initialization, optimization
%   will be initialized using the value of unconstrained parameters
%   contained in this object. The initial value of sigma in Psi will be
%   ignored. When optimization completes, Psi will contain the optimal
%   theta and sigma.
        Psi
        
%FitMethod - The method used to fit the LME model - either 'ML' for maximum 
%   likelihood or 'REML' for restricted maximum likelihood.
        FitMethod                    
    end

% Public properties storing size information.    
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
  
% Public properties storing optimization info. 
    properties(GetAccess=public, SetAccess=protected)

%Optimizer - A string containing the name of the optimizer.        
        Optimizer = 'quasinewton';
        
%OptimizerOptions - A structure containing the optimizer options.        
        OptimizerOptions = struct([]);
        
%CheckHessian - A logical scalar indicating whether to perform the Hessian 
%   checks in initstats.    
        CheckHessian = false;
        
    end
    
% Public property storing the initialization method for optimization.
    properties(GetAccess=public, SetAccess=public)

%InitializationMethod - A string containing the method used to compute the 
%   initial value of theta to start the optimization.        
        InitializationMethod = 'default';
        
    end
    
% Public properties for estimated LME model parameters.   
    properties (GetAccess=public, SetAccess=protected)
%betaHat - Estimated fixed effects vector.
        betaHat
        
%bHat - Estimated Best Linear Unbiased Predictor (EBLUP) of random effects 
%   vector.
        bHat
        
%DeltabHat - Estimated normalized posterior mode of random effects vector.
        DeltabHat        
        
%sigmaHat - Estimated residual standard deviation.
        sigmaHat
        
%thetaHat - Estimated vector of unconstrained parameters for matrix D.
        thetaHat        
        
%loglikHat - Maximized log likelihood (if FitMethod = 'ML') or maximized 
%   restricted log likelihood (if FitMethod = 'REML').        
        loglikHat                    
    end
        
% Public properties for stats info on LME model parameters.
    properties (GetAccess=public, SetAccess=protected)        
%covbetaHat - Estimated covariance of betaHat.        
        covbetaHat
                
%covthetaHatlogsigmaHat - Estimated covariance of [thetaHat;log(sigmaHat)].        
        covthetaHatlogsigmaHat      
        
%covetaHatlogsigmaHat - Estimated covariance of [etaHat;log(sigmaHat)]. 
%   etaHat is the Natural parameter vector corresponding to thetaHat.        
        covetaHatlogsigmaHat
        
%covbetaHatbHat - Estimated covariance of [betaHat-beta;bHat-b].
        covbetaHatbHat
    end        
    
% Public properties storing whether sigma is fixed and if so its value.
    properties (GetAccess=public, SetAccess=protected)
%isSigmaFixed - True if residual standard deviation sigma is fixed.
        isSigmaFixed = false;
        
%sigmaFixed - Scalar fixed value for residual standard deviation.
        sigmaFixed = NaN;
    end    
    
% Public properties storing state information.
    properties (GetAccess=public, SetAccess=protected)
%isFitToData - True if this is a fitted StandardLinearMixedModel object.
        isFitToData = false;
        
%isReadyForStats - True if this StandardLinearMixedModel object is ready 
%   for stats.
        isReadyForStats = false;
    end
     
% Internal constants.
    properties (Constant=true, Hidden=true)       
        
%AllowedFitMethods - The fitting methods that we currently support.        
        AllowedFitMethods = {'ML','REML'};                        
        
    end        
  
% Private properties storing precomputed matrix/vector products.    
    properties (Access=private)
%XtX - Precomputed value of X'*X where X is the fixed effects design
%   matrix.
        XtX
%Xty - Precomputed value of X'*y where X is the fixed effects design
%   matrix and y is the response vector.
        Xty
        
%XtZ - Precomputed value of X'*Z where X is the fixed effects design matrix
%   and Z is the random effects design matrix.
        XtZ
        
%Zty - Precomputed value of Z'*y where Z is the random effects design
%   matrix and y is the response vector.
        Zty

%ZtZ - Precomputed value of Z'*Z where Z is the random effects design
%   matrix.
        ZtZ        
    end
    
% Set methods for X, y, Z, Psi, FitMethod and InitializationMethod.
    methods
        
        function slme = set.y(slme,newy)
        
            if ~isempty(slme.y)
                % (1) Validate newy.
                newy = validatey(slme,newy);
                
                % (2) Invalidate the fit.
                slme = invalidateFit(slme);
            end
            
            % (3) Set the new y.
            slme.y = newy;
            
        end % end of set.y
        
        function slme = set.X(slme,newX)
            
            if ~isempty(slme.X)
                % (1) Validate newX.
                newX = validateX(slme,newX);
                
                % (2) Invalidate the fit.
                slme = invalidateFit(slme);
            end
            
            % (3) Set the new X.
            slme.X = newX;
            
        end % end of set.X
        
        function slme = set.Z(slme,newZ)
            
            if ~isempty(slme.Z)
                % (1) Validate newZ.
                newZ = validateZ(slme,newZ);
                
                % (2) Invalidate the fit.
                slme = invalidateFit(slme);
            end            
            
            % (3) Set the new Z.
            slme.Z = newZ;
            
        end % end of set.Z
        
        function slme = set.Psi(slme,newPsi)
            
            if ~isempty(slme.Psi)
                % (1) Validate newPsi.
                newPsi = validatePsi(slme,newPsi);
                
                % (2) Invalidate the fit.
                slme = invalidateFit(slme);
            end
            
            % (3) Set the new Psi.
            slme.Psi = newPsi;
            
        end % end of set.Psi
        
        function slme = set.FitMethod(slme,newFitMethod)
            
            if ~isempty(slme.FitMethod)
                % (1) Validate newFitMethod.
                newFitMethod = validateFitMethod(slme,newFitMethod);
                
                % (2) Invalidate the fit.
                slme = invalidateFit(slme);
            end
            
            % (3) Set the new FitMethod.
            slme.FitMethod = newFitMethod;            
            
        end % end of set.FitMethod
        
        function slme = set.InitializationMethod(slme,newInitializationMethod)
            
            if ~isempty(slme.InitializationMethod)
                % (1) Validate newInitializationMethod.
                newInitializationMethod = validateInitializationMethod(slme,newInitializationMethod);
                
                % (2) Invalidate the fit.
                slme = invalidateFit(slme);
            end
            
            % (3) Set the new InitializationMethod.
            slme.InitializationMethod = newInitializationMethod;
            
        end % end of set.InitializationMethod.
        
    end
    
% Protected methods for input validation.
    methods (Access=protected) 
        
        function FitMethod = validateFitMethod(slme,FitMethod) %#ok<INUSL>
            
            % (1) FitMethod is a string containing either 'ML' or 'REML'.
            FitMethod = internal.stats.getParamVal(FitMethod,classreg.regr.lmeutils.StandardLinearMixedModel.AllowedFitMethods,'FitMethod');
            
        end % end of validateFitMethod.        
        
    end
    
% Protected methods related to fitting.
    methods (Access=protected)        
        
        function fun = makeObjectiveFunctionForMinimization(slme)
%makeObjectiveFunctionForMinimization - Get objective function to minimize.
%   fun = makeObjectiveFunctionForMinimization(slme) takes an object slme
%   of type StandardLinearMixedModel and depending on slme.FitMethod makes
%   either the objective function for ML or REML.

            % (1) Make an objective function for minimization. 
            %     (a) If slme.FitMethod is 'ML' this is the negative 
            %     profiled log likelihood (beta and sigma profiled out).
            %
            %     (b) If slme.FitMethod is 'REML' this is the negative 
            %     profiled restricted log likelihood (sigma profiled out).
            switch lower(slme.FitMethod)
                case 'ml'
                    fun = makeNegBetaSigmaProfiledLogLikelihoodAsFunctionOfTheta(slme);
                case 'reml'
                    fun = makeNegSigmaProfiledRestrictedLogLikelihoodAsFunctionOfTheta(slme);
                otherwise 
                    % <entry key="BadFitMethod">''FitMethod'' must be either ''ML'' or ''REML''.</entry>
                    error(message('stats:classreg:regr:lmeutils:StandardLinearMixedModel:BadFitMethod'));
            end
            
        end % makeObjectiveFunctionForMinimization.
        
        function theta0 = initializeTheta(slme)
%initializeTheta - Return an initial value of theta to start optimization.            
            
            % If initialization method is:
            %   (a) 'default' - then extract theta0 from slme.Psi
            %   (b) 'random'  - then generate theta0 from a standard Normal            
            switch lower(slme.InitializationMethod)
                case 'default'
                    % Use the value stored in slme.Psi.
                    theta0 = getUnconstrainedParameters(slme.Psi);
                case 'random'
                    % Use a random value from N(0,1).
                    theta0 = getUnconstrainedParameters(slme.Psi);
                    theta0 = randn(length(theta0),1);
                otherwise
                    % <entry key="BadInitializationMethod">''InitializationMethod'' must be either ''Default'' or ''Random''.</entry>
                    error(message('stats:classreg:regr:lmeutils:StandardLinearMixedModel:BadInitializationMethod'));
            end
                        
        end % end of initializeTheta.
        
        function warnAboutPerfectFit(slme,theta0)
%warnAboutPerfectFit - Warn if theta0 gives a nearly perfect fit.            
%   warnAboutPerfectFit(slme,theta0) takes a StandardLinearMixedModel 
%   object slme and an initial value of the unconstrained parameter vector
%   theta0 and warns the user if we have a nearly perfect fit using theta0.
%   We simply solve the mixed model equations at theta0 and if the
%   estimated residual standard deviation is "small", we display a warning
%   message.

            % (1) Solve mixed model equations at theta0.
            [~,~,sigmaHat,~,~,~] = solveMixedModelEquations(slme,theta0);

            % (2) Warn about perfect fit if sigmaHat is "small".
            sigmaHatTol = sqrt(eps(class(sigmaHat)));
            if abs(sigmaHat) <= sigmaHatTol*std(slme.y)
                % <entry key="Message_PerfectFit">Initial estimate of residual noise standard deviation is almost zero. This may indicate a nearly perfect fit to data.</entry>
                warning(message('stats:classreg:regr:lmeutils:StandardLinearMixedModel:Message_PerfectFit'));
            end
            
        end % end of warnAboutPerfectFit.
                
        function thetaHat = solveForThetaHat(slme)
%solveForThetaHat - Find the optimal value of theta and return it.
%   thetaHat = solveForThetaHat(slme) takes a StandardLinearMixedModel 
%   object slme, solves for the optimal theta and returns it. theta is the
%   unconstrained vector that parameterizes the matrix D.

            % (1) Make an objective function for minimization. 
            fun = makeObjectiveFunctionForMinimization(slme);                      
            
            % (2) Get theta0, the initial value of theta.
            theta0 = initializeTheta(slme);
                    
            % (3) Warn about perfect fit at theta0.
            warnAboutPerfectFit(slme,theta0);
            
            % (4) Minimize fun, starting at theta0. Use the optimization
            % routine and optimization options stored in object slme.            
            [thetaHat,cause] = doMinimization(slme,fun,theta0);
            
            % (5) If cause = 0 or 1 then we are good, otherwise display a
            % message that we did not converge.
            if ( cause ~= 0 && cause ~= 1 )
                % <entry key="Message_OptimizerUnableToConverge">Optimizer {0} was unable to converge to a solution.</entry>
                warning(message('stats:classreg:regr:lmeutils:StandardLinearMixedModel:Message_OptimizerUnableToConverge',slme.Optimizer));                         
            end
            
        end % end of solveForThetaHat.                                        
        
        function H = objectiveFunctionHessianAtThetaHat(slme)
%objectiveFunctionHessianAtThetaHat - Objective function Hessian.           
%   H = objectiveFunctionHessianAtThetaHat(slme) takes a fitted object slme
%   of type StandardLinearMixedModel and gets the Hessian of the objective
%   function used for minimization evaluated at the optimal slme.thetaHat.

            % (1) Make an objective function for minimization.
            fun = makeObjectiveFunctionForMinimization(slme);

            % (2) Get the Hessian of fun at slme.thetaHat.
            wantRegularized = false;
            H = slme.getHessian(fun,slme.thetaHat,wantRegularized);         

        end % end of objectiveFunctionHessianAtThetaHat.
        
        function checkObjectiveFunctionHessianAtThetaHat(slme)
%checkObjectiveFunctionHessianAtThetaHat - Positive definiteness check.
%   checkObjectiveFunctionHessianAtThetaHat(slme) takes a fitted object
%   slme of type StandardLinearMixedModel and checks if the Hessian H of 
%   the objective function used for minimization evaluated at the optimal 
%   slme.thetaHat is positive definite or not. If H is positive definite 
%   then all is well, otherwise a warning message is printed to the screen.

            % (1) Get Hessian of objective function for minimization at
            % slme.thetaHat.
            H = objectiveFunctionHessianAtThetaHat(slme);
            
            % (2) Get the warning message string if H is not positive
            % definite.
            switch lower(slme.FitMethod)                    
                case 'ml'
                    % <entry key="Message_NotSPDHessian_ML">The Hessian of negative profiled log likelihood as a function of unconstrained parameters is not positive definite at convergence.</entry>
                    msgID = 'stats:classreg:regr:lmeutils:StandardLinearMixedModel:Message_NotSPDHessian_ML';
                case 'reml'
                    % <entry key="Message_NotSPDHessian_REML">The Hessian of negative profiled restricted log likelihood as a function of unconstrained parameters is not positive definite at convergence.</entry>   
                    msgID = 'stats:classreg:regr:lmeutils:StandardLinearMixedModel:Message_NotSPDHessian_REML';
                otherwise
                    % <entry key="BadFitMethod">''FitMethod'' must be either ''ML'' or ''REML''.</entry>      
                    error(message('stats:classreg:regr:lmeutils:StandardLinearMixedModel:BadFitMethod'));
            end
            
            % (3) Ensure H is positive definite. If not, display the
            % message msg created above.            
            slme.checkPositiveDefinite(H,msgID);
            
        end % end of checkObjectiveFunctionHessianAtThetaHat.        
        
        function [betaHat,bHat,sigmaHat,r2,R,R1,Deltab] = solveMixedModelEquations(slme,theta)
%solveMixedModelEquations - Solve mixed model equations.
%   [betaHat,bHat,sigmaHat,r2,R,R1] = solveMixedModelEquations(slme,theta)
%   takes a StandardLinearMixedModel object slme and solves the mixed model
%   equations at the specified value of theta. First, we solve mixed model
%   equations using theta to get betaHat and bHat. Then we get sigmaHat
%   (for ML or REML), r2, R and R1. The meaning of r2, R and R1 is defined
%   in internal LME specs. Deltab is the estimated normalized posterior
%   mode of the random effects vector.
    
            % (1) Get X, y, Z, N, p and q from slme. Also get various
            % matrix and vector products to reuse precomputed values.
            X = slme.X;
            y = slme.y;
            Z = slme.Z;
            N = slme.N;
            p = slme.p;
            q = slme.q;            
            XtX = slme.XtX;
            Xty = slme.Xty;
            XtZ = slme.XtZ;
            Zty = slme.Zty;
            ZtZ = slme.ZtZ;
            
            % (2) Get slme.Psi, set current theta and set sigma = 1.
            Psi = slme.Psi;
            Psi = setUnconstrainedParameters(Psi,theta);
            Psi = setSigma(Psi,1);

            % (3) Get lower triangular Cholesky factor of D matrix in Psi.
            % Lambda has size q by q where q = size(Z,2).
            Lambda = getLowerTriangularCholeskyFactor(Psi);

            % (4) Compute the U matrix. U is N by q.
            %U = Z*Lambda;
            
            % TODO: Try to make steps 5 onward more efficient.
            % (5) Get R and S such that S*R'*R*S' = U'*U + eye(q) where R
            % and S are q by q. S is a permutation matrix: S*S' = eye(q).
            Iq = spdiags(ones(q,1),0,q,q);
            %[R,status,S] = chol(U'*U + Iq);
            [R,status,S] = chol(sparse(Lambda'*ZtZ*Lambda + Iq));
            
            % (6) Make sure status is zero.
            if (status ~= 0)
                % <entry key="ErrorSparseCholesky">An error occured during sparse Cholesky factorization.</entry>
                error(message('stats:classreg:regr:lmeutils:StandardLinearMixedModel:ErrorSparseCholesky'));
            end
            
            % For steps (7) onwards, we may get a rank deficient X or a
            % singular factor R1. Hence we disable the nearly singular
            % matrix warning and then enable it at the end of this method
            % or on error.
            warnState = warning('query','all');
            warning('off','MATLAB:nearlySingularMatrix');
            warning('off','MATLAB:illConditionedMatrix');
            warning('off','MATLAB:singularMatrix');
            cleanupObj = onCleanup(@() warning(warnState));
            
            % (7) Compute P1 (q by q), Q1 (p by q) and R1 (p by p).
            %P1 = S*R';
            %Q1 = ((X'*U)*S) / R;
            Q1 = ((XtZ*Lambda)*S) / R;
            
            R1R1t = XtX - Q1*Q1';
            try
                %R1 = chol(X'*X - Q1*Q1','lower');
                %R1 = chol(XtX - Q1*Q1','lower');
                R1 = chol(R1R1t,'lower');
            catch ME %#ok<NASGU>
                % We know from theory that R1 must exist.
                %R1 = classreg.regr.lmeutils.covmats.CovarianceMatrix.singularLowerChol(X'*X - Q1*Q1');
                R1 = classreg.regr.lmeutils.covmats.CovarianceMatrix.singularLowerChol(R1R1t);
            end
                        
            % (8) Compute cDeltab (q by 1) and cbeta (p by 1).
            %cDeltab = P1 \ (U'*y);
            %cDeltab = R' \ (S'*(U'*y));
            cDeltab = R' \ (S'*((Lambda'*Zty)));
            %cbeta = R1 \ (X'*y - Q1*cDeltab);
            cbeta = R1 \ (Xty - Q1*cDeltab);
            
            % (9) Compute betaHat and Deltab.
            betaHat = R1' \ cbeta;
            %Deltab = P1' \ (cDeltab - Q1'*betaHat);
            Deltab = S*(R \ (cDeltab - Q1'*betaHat));

            % (10) Compute bHat.
            bHat = Lambda * Deltab;
            
            % (11) Compute r2.
            r2 = sum(Deltab.^2) + sum((y - X*betaHat - Z*bHat).^2);
            
            % (12) Compute sigmaHat.
            if slme.isSigmaFixed
                sigmaHat = slme.sigmaFixed;
            else
                switch lower(slme.FitMethod)
                    case 'ml'
                        sigma2 = r2/N;
                        sigmaHat = sqrt(sigma2);
                    case 'reml'
                        sigma2 = r2/(N - p);
                        sigmaHat = sqrt(sigma2);
                    otherwise
                        % <entry key="BadFitMethod">''FitMethod'' must be either ''ML'' or ''REML''.</entry>
                        error(message('stats:classreg:regr:lmeutils:StandardLinearMixedModel:BadFitMethod'));
                end
            end
            % NOTE: Restore the warning state to that saved in warnState.
            % This is automatically done by the cleanupObj going out of
            % scope. warning(warnState);
            
        end % solveMixedModelEquations.           
    
        function [pred,varpred] = getEstimateAndVariance(slme,Xnew,Znew,betaBar,DeltabBar,theta,sigma,type)
%[pred,varpred] = getEstimateAndVariance(slme,Xnew,Znew,betaBar,DeltabBar,theta,sigma,type)
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
            assert( p == slme.p && q == slme.q );
            assert( size(betaBar,1) == p && size(betaBar,2) == 1 );
            assert( size(DeltabBar,1) == q && size(DeltabBar,2) == 1 );
            
            % 2. Results should be insensitive to sigma if it is fixed.
            if slme.isSigmaFixed
                sigma = slme.sigmaFixed;
            end

            % 3. Get Lambda.
            Psi = slme.Psi;
            Psi = Psi.setUnconstrainedParameters(theta);
            Lambda = getLowerTriangularCholeskyFactor(Psi);
                
            % 4. Get pred.                
            pred = Xnew*betaBar + Znew*(Lambda*DeltabBar);
            pred = full(pred);
            
            % 5. Need varpred?
            if nargout > 1
                warnState = warning('query','all');
                warning('off','MATLAB:nearlySingularMatrix');
                warning('off','MATLAB:illConditionedMatrix');
                warning('off','MATLAB:singularMatrix');
                cleanupObj = onCleanup(@() warning(warnState));
                
                % 5.1 Compute R, S, Q1 and R1. For LMEs, we know for sure
                % that R1 must exist.
                X = slme.X;
                Z = slme.Z;
                U = Z*Lambda;
                Iq = spdiags(ones(q,1),0,q,q);
                [R,status,S] = chol(U'*U + Iq);
                if ( status ~= 0 )
                    error(message('stats:classreg:regr:lmeutils:StandardLinearMixedModel:ErrorSparseCholesky'));
                end
                Q1 = (X'*U*S) / R;
                try
                    R1 = chol(X'*X - Q1*Q1','lower');
                catch ME %#ok<NASGU>
                    % We know from theory that R1 must exist.
                    R1 = classreg.regr.lmeutils.covmats.CovarianceMatrix.singularLowerChol(X'*X - Q1*Q1');
                end
                
                % 5.2 Compute Cb and Ca.
                if ( M <= q )
                    Cb = R' \ (S'*(Lambda'*Znew'));
                else
                    Cb = (R' \ (S'*Lambda')) * Znew';
                end
                Ca = R1 \ (Xnew' - Q1*Cb);
                
                % 5.3 Compute varpred and make it into a column vector.
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

% Protected methods to "make" the various likelihood functions.
    methods (Access=protected)  
        
        function fun = makeNegBetaSigmaProfiledLogLikelihoodAsFunctionOfTheta(slme)
%makeNegBetaSigmaProfiledLogLikelihoodAsFunctionOfTheta - This returns the
%   negative profiled log likelihood with beta and sigma profiled out as a
%   function of theta.

            fun = @f0;            
            function y0 = f0(theta)                
                 L = BetaSigmaProfiledLogLikelihood(slme,theta);
                y0 = -1*L;
                % Don't let y0 be equal to -Inf.
                y0 = max(-realmax,y0);
            end
            
        end % end of makeNegBetaSigmaProfiledLogLikelihoodAsFunctionOfTheta.
       
        function fun = makeNegSigmaProfiledRestrictedLogLikelihoodAsFunctionOfTheta(slme)
%makeNegSigmaProfiledRestrictedLogLikelihoodAsFunctionOfTheta - This
%   returns the negative profiled restricted log likelihood with sigma
%   profiled out as a function of theta.

            fun = @f1;
            function y1 = f1(theta)
                 L = SigmaProfiledRestrictedLogLikelihood(slme,theta);
                y1 = -1*L;
                % Don't let y1 be equal to -Inf.
                y1 = max(-realmax,y1);
            end
            
        end % end of makeNegSigmaProfiledRestrictedLogLikelihoodAsFunctionOfTheta.
        
        function fun = makeNegBetaProfiledLogLikelihoodAsFunctionOfThetaLogSigma(slme)
%makeNegBetaProfiledLogLikelihoodAsFunctionOfThetaLogSigma - This returns
%   the negative profiled log likelihood with beta profiled out as a
%   function of x = [theta;log(sigma)].
           
            fun = @f2;
            function y2 = f2(x)
                % (1) Extract theta and logsigma from x.
                   theta = x(1:end-1);
                logsigma = x(end);
                % (2) Evaluate BetaProfiledLogLikelihood at (theta,sigma)
                % and negate the result.
                   sigma = exp(logsigma);
                 L = BetaProfiledLogLikelihood(slme,theta,sigma);
                y2 = -1*L;
            end
            
        end % end of makeNegBetaProfiledLogLikelihoodAsFunctionOfThetaLogSigma.
        
        function fun = makeNegRestrictedLogLikelihoodAsFunctionOfThetaLogSigma(slme)
%makeNegRestrictedLogLikelihoodAsFunctionOfThetaLogSigma - This returns the
%   negative restricted log likelihood as a function of x =
%   [theta;log(sigma)].
            
            fun = @f3;
            function y3 = f3(x)
                % (1) Extract theta and logsigma from x.
                   theta = x(1:end-1);
                logsigma = x(end);
                % (2) Evaluate RestrictedLogLikelihood at (theta,sigma) and
                % negate the result.
                   sigma = exp(logsigma);
                   L = RestrictedLogLikelihood(slme,theta,sigma);
                  y3 = -1*L;                
            end
            
        end % end of makeNegRestrictedLogLikelihoodAsFunctionOfThetaLogSigma.
        
        function fun = makeNegBetaProfiledLogLikelihoodAsFunctionOfEtaLogSigma(slme)
%makeNegBetaProfiledLogLikelihoodAsFunctionOfEtaLogSigma - This returns
%   the negative profiled log likelihood with beta profiled out as a
%   function of x = [eta;log(sigma)] where eta is the Natural parameter
%   vector.

            fun = @f6;
            function y6 = f6(x)
                % (1) Extract eta, logsigma and sigma from x.
                     eta = x(1:end-1);
                logsigma = x(end);
                if slme.isSigmaFixed
                    sigma = slme.sigmaFixed;
                else
                    sigma = exp(logsigma);
                end
                % (2) Set sigma and eta in slme.Psi.
                Psi = slme.Psi;
                Psi = setSigma(Psi,sigma);
                Psi = setNaturalParameters(Psi,eta);
                % (3) Get theta from Psi.
                theta = getUnconstrainedParameters(Psi);
                % (4) Evaluate BetaProfiledLogLikelihood at (theta,sigma)
                % and negate the result.
                L = BetaProfiledLogLikelihood(slme,theta,sigma);
                y6 = -1*L;                
            end

        end % end of makeNegBetaProfiledLogLikelihoodAsFunctionOfEtaLogSigma.
        
        function fun = makeNegRestrictedLogLikelihoodAsFunctionOfEtaLogSigma(slme)
%makeNegRestrictedLogLikelihoodAsFunctionOfEtaLogSigma - This returns the
%   negative restricted log likelihood as a function of x =
%   [eta;log(sigma)] where eta is the Natural parameter vector.
            
            fun = @f7;
            function y7 = f7(x)
                % (1) Extract eta, logsigma and sigma from x.
                     eta = x(1:end-1);
                logsigma = x(end);
                if slme.isSigmaFixed
                    sigma = slme.sigmaFixed;
                else
                    sigma = exp(logsigma);
                end
                % (2) Set sigma and eta in slme.Psi.
                Psi = slme.Psi;
                Psi = setSigma(Psi,sigma);
                Psi = setNaturalParameters(Psi,eta);
                % (3) Get theta from Psi.
                theta = getUnconstrainedParameters(Psi);
                % (4) Evaluate RestrictedLogLikelihood at (theta,sigma) and
                % negate the result.
                L = RestrictedLogLikelihood(slme,theta,sigma);
                y7 = -1*L;
            end

            
        end % end of makeNegRestrictedLogLikelihoodAsFunctionOfEtaLogSigma.
        
    end
    
% Protected methods related to computing various likelihoods. 
    methods (Access=protected)    
        
        function LogLik = LogLikelihood(slme,theta,sigma,beta)
%LogLikelihood - Compute the log likelihood of a StandardLinearMixedModel.
%   LogLik = LogLikelihood(slme,theta,sigma,beta) returns the log
%   likelihood of the StandardLinearMixedModel slme at the specified values
%   of theta, sigma and beta.

            if slme.isSigmaFixed
                sigma = slme.sigmaFixed;
            end

            % (1) Get X, y, Z, N and q from slme.
            X = slme.X;
            y = slme.y;
            Z = slme.Z;
            N = slme.N;
            q = slme.q;
           
            % (2) Get slme.Psi, set current theta and set sigma = 1.
            Psi = slme.Psi;
            Psi = setUnconstrainedParameters(Psi,theta);
            Psi = setSigma(Psi,1);

            % (3) Get lower triangular Cholesky factor of D matrix in Psi.
            % Lambda has size q by q where q = size(Z,2).
            Lambda = getLowerTriangularCholeskyFactor(Psi);

            % (4) Compute the U matrix. U is N by q.
            U = Z*Lambda;
            
            % TODO: Try to make steps 5 onward more efficient.
            % (5) Get R and S such that S*R'*R*S' = U'*U + eye(q) where R
            % and S are q by q. S is a permutation matrix: S*S' = eye(q).
            Iq = spdiags(ones(q,1),0,q,q);
            [R,status,S] = chol(U'*U + Iq);
 
            % (6) Make sure status is zero.
            if (status ~= 0)
                % <entry key="ErrorSparseCholesky">An error occured during sparse Cholesky factorization.</entry>
                error(message('stats:classreg:regr:lmeutils:StandardLinearMixedModel:ErrorSparseCholesky'));
            end            
            
            % (7) Compute P1 (q by q).
            P1 = S*R';

            % (8) Compute Deltab for fixed beta.
            Deltab = P1' \ (P1 \ (U'*(y-X*beta)) );
            
            % (9) Compute b.
            b = Lambda * Deltab;
            
            % (10) Compute r2.
            r2 = sum(Deltab.^2) + sum((y - X*beta - Z*b).^2);
            
            % (11) Get log abs determinant of upper triangular matrix R.
            logAbsDetR = slme.logAbsDetTriangular(R);

            % (12) Compute the log likelihood.
            sigma2 = sigma^2;
            LogLik = (-N/2)*log(2*pi*sigma2) - r2/(2*sigma2) - logAbsDetR;            

        end % end of LogLikelihood.
        
        function PLogLik = BetaSigmaProfiledLogLikelihood(slme,theta)
%BetaSigmaProfiledLogLikelihood - Compute the profiled log likelihood with beta and sigma profiled out.
%   PLogLik = BetaSigmaProfiledLogLikelihood(slme,theta) returns the 
%   profiled log likelihood of a StandardLinearMixedModel slme at the 
%   specified value of theta. Both beta and sigma have been profiled out.           

            if slme.isSigmaFixed
                PLogLik = BetaProfiledLogLikelihood(slme,theta,slme.sigmaFixed);
            else
                % (1) Using the given value of theta, solve mixed model 
                % equations to get r2 and R.
                [~,~,~,r2,R,~] = solveMixedModelEquations(slme,theta);

                % (2) Get log abs determinant of upper triangular matrix R.
                logAbsDetR = slme.logAbsDetTriangular(R);

                % (3) Compute profiled log likelihood.
                N = slme.N;
                PLogLik = (-N/2)*( 1 + log( 2*pi*r2/N ) ) - logAbsDetR;            
            end
            
        end % end of BetaSigmaProfiledLogLikelihood.
        
        function PLogLik = BetaProfiledLogLikelihood(slme,theta,sigma)
%BetaProfiledLogLikelihood - Compute the profiled log likelihood with beta profiled out.
%   PLogLik = BetaProfiledLogLikelihood(slme,theta,sigma) returns the
%   profiled log likelihood of a StandardLinearMixedModel slme at the
%   specified values of theta and sigma with beta profiled out.

            if slme.isSigmaFixed
                sigma = slme.sigmaFixed;
            end

            % (1) Using the given value of theta, solve mixed model 
            % equations to get r2 and R.
            [~,~,~,r2,R,~] = solveMixedModelEquations(slme,theta);
            
            % (2) Get log abs determinant of upper triangular matrix R.
            logAbsDetR = slme.logAbsDetTriangular(R);

            % (3) Compute profiled log likelihood (sigma not profiled out).
            N = slme.N;
            sigma2 = sigma^2;
            PLogLik = (-N/2)*log(2*pi*sigma2) - r2/(2*sigma2) - logAbsDetR;
            
        end % end of BetaProfiledLogLikelihood.                
        
        function PRLogLik = SigmaProfiledRestrictedLogLikelihood(slme,theta)
%SigmaProfiledRestrictedLogLikelihood - Compute the profiled restricted log likelihood.
%   PRLogLik = SigmaProfiledRestrictedLogLikelihood(slme,theta) returns the
%   profiled restricted log likelihood of a StandardLinearMixedModel slme
%   at the specified value of theta with sigma profiled out.

            if slme.isSigmaFixed
                PRLogLik = RestrictedLogLikelihood(slme,theta,slme.sigmaFixed);
            else
                % (1) Using the given value of theta, solve mixed model
                % equations to get r2, R and R1.
                [~,~,~,r2,R,R1] = solveMixedModelEquations(slme,theta);

                % (2) Get log abs determinant of upper triangular R and lower
                % triangular R1.
                 logAbsDetR = slme.logAbsDetTriangular(R);
                logAbsDetR1 = slme.logAbsDetTriangular(R1);

                % (3) Compute profiled restricted log likelihood.
                N = slme.N;
                p = slme.p;
                PRLogLik = (-(N-p)/2) * ( 1 + log( 2*pi*r2/(N-p) ) ) ...
                    - logAbsDetR - logAbsDetR1;            
            end
        end % end of SigmaProfiledRestrictedLogLikelihood.        
        
        function RLogLik = RestrictedLogLikelihood(slme,theta,sigma)
%RestrictedLogLikelihood - Compute the restricted log likelihood. 
%   RLogLik = RestrictedLogLikelihood(slme,theta,sigma) returns the
%   restricted log likelihood of a StandardLinearMixedModel at the
%   specified values of theta and sigma.

            if slme.isSigmaFixed
                sigma = slme.sigmaFixed;
            end

            % (1) Using the given value of theta, solve mixed model
            % equations to get r2, R and R1.
            [~,~,~,r2,R,R1] = solveMixedModelEquations(slme,theta);
            
            % (2) Get log abs determinant of upper triangular R and lower
            % triangular R1.
             logAbsDetR = slme.logAbsDetTriangular(R);
            logAbsDetR1 = slme.logAbsDetTriangular(R1);

            % (3) Compute restricted log likelihood.
            N = slme.N;
            p = slme.p;
            sigma2 = sigma^2;
            RLogLik = (-(N-p)/2)*log(2*pi*sigma2) - r2/(2*sigma2) ...
                - logAbsDetR - logAbsDetR1;
            
        end % end of RestrictedLogLikelihood.                                                              
        
    end       
        
% Covariance of betaHat and [betaHat - beta;bHat - b] as functions of theta and sigma. 
    methods (Access=protected) 
        
        function covbetaHat = covBetaHatAsFunctionOfThetaSigma(slme,theta,sigma)
%covBetaHatAsFunctionOfThetaSigma - Covariance matrix of betaHat as function of theta and sigma.
%   covbetaHat = covBetaHatAsFunctionOfThetaSigma(slme,theta,sigma) returns
%   the covariance matrix of betaHat as a function of theta and sigma.
            
            if slme.isSigmaFixed
                sigma = slme.sigmaFixed;
            end
            
            % (1) Given theta, get R1.
            [~,~,~,~,~,R1] = solveMixedModelEquations(slme,theta);
            
            % (2) Return covbetaHat = sigma^2 * inv(R1*R1') 
            %                       = sigma^2 * inv(R1)'*inv(R1)
            sigma2 = sigma^2;
            invR1 = R1 \ eye(size(R1));
            covbetaHat = sigma2*(invR1'*invR1);
            
        end % end of covBetaHatAsFunctionOfThetaSigma.
        
        function f = makecovcTBetaHatAsFunctionOfThetaLogSigma(slme,c)   
%   f = makecovcTBetaHatAsFunctionOfThetaLogSigma(slme,c) takes a p by 1
%   vector c and gets the covariance of c'*betaHat as a function of
%   [theta;log(sigma)].

                f = @f4;
                function ret = f4(x)
                    theta = x(1:end-1);
                    logsigma = x(end);
                    sigma = exp(logsigma);                    
                    ret = c'*covBetaHatAsFunctionOfThetaSigma(slme,theta,sigma)*c;                    
                end                
                
        end % makecovcTBetaHatAsFunctionOfThetaLogSigma.
                       
        function covbetaHatbHat = covBetaHatBHatAsFunctionOfThetaSigma(slme,theta,sigma)
%   covbetaHatbHat = covBetaHatBHatAsFunctionOfThetaSigma(slme,theta,sigma) 
%   returns the covariance of [betaHat - beta;bHat - b] as a function of
%   theta and sigma.

            if slme.isSigmaFixed
                sigma = slme.sigmaFixed;
            end

            % (1) Get X, Z, p and q
            X = slme.X;
            Z = slme.Z;
            p = slme.p;
            q = slme.q;
            
            % (2) Get Psi from slme and set its unconstrained parameters.
            Psi = slme.Psi;
            Psi = Psi.setUnconstrainedParameters(theta);
            
            % (3) Get lower triangular Cholesky factor from Psi.
            Lambda = getLowerTriangularCholeskyFactor(Psi);
            
            % (4) Get U.
            U = Z*Lambda;
            
            % (5) Form the matrix [X'*X X'*U;U'*X U'*U + I].
            M = zeros(p+q);
            M(1:p,1:p) = X'*X;
            M(1:p,p+1:end) = X'*U;
            M(p+1:end,1:p) = M(1:p,p+1:end)';
            M(p+1:end,p+1:end) = U'*U + eye(q);
            
            % (6) Get Cholesky factor C such that C*C' = M.
            try 
                C = chol(M,'lower');
            catch ME %#ok<NASGU>
                C = classreg.regr.lmeutils.covmats.CovarianceMatrix.singularLowerChol(M);
            end
            
            % (7) Get G = [I_p 0;0 Lambda^T].
            G = zeros(p+q);
            G(1:p,1:p) = eye(p);
            G(p+1:end,p+1:end) = Lambda';
                        
            % (8) Get T = inv(C) * G.
            T = C \ G;
            
            % (9) Get sigma^2 * (T'*T) as the required covariance.
            covbetaHatbHat = (sigma^2)*(T'*T);
                    
        end % end of covBetaHatBHatAsFunctionOfThetaSigma.        
        
        function f = makecovcTBetaHatBHatAsFunctionOfThetaLogSigma(slme,c)
%   f = makecovcTBetaHatBHatAsFunctionOfThetaLogSigma(slme,c) takes a (p+q)
%   by 1 vector c and gets the covariance of c'*[betaHat;bHat] as a
%   function of [theta;log(sigma)].

                f = @f5;
                function ret = f5(x)
                    theta = x(1:end-1);
                    logsigma = x(end);
                    sigma = exp(logsigma);
                    %ret = c'*covBetaHatBHatAsFunctionOfThetaSigma(slme,theta,sigma)*c;
                    
                    % For a LME model, the joint covariance of betaHat -
                    % beta and bHat - b does depends only on theta and
                    % sigma. Since we don't really want pred from 
                    % getEstimateAndVariance, we can set betaBat and
                    % DeltabBat equal to NaNs.
                    p = slme.p;
                    q = slme.q;
                    Xnew = c(1:p,1)';
                    Znew = c(p+1:p+q,1)';
                    betaBar = NaN(p,1);
                    DeltabBar = NaN(q,1);
                    [~,ret] = getEstimateAndVariance(slme,Xnew,Znew,betaBar,DeltabBar,theta,sigma);                    
                end
            
        end % end of makecovcTBetaHatBHatAsFunctionOfThetaLogSigma.            
        
        % Get estimated covariance of [betaHat - beta;bHat - b].
        function C = covBetaHatBHat(slme)
%covBetaHatBHat - Estimated covariance matrix of [betaHat-beta; bHat-b].          
%   C = covBetaHatBHat(slme) returns the estimated covariance matrix of
%   [betaHat - beta; bHat - b].
            
            % Calculate the covariance if not already stored in object.
            if isempty(slme.covbetaHatbHat)
                C = covBetaHatBHatAsFunctionOfThetaSigma(slme,...
                    slme.thetaHat,slme.sigmaHat);
            else
                C = slme.covbetaHatbHat;
            end            
            
        end % end of covBetaHatBHat.
        
    end
           
% Protected methods to get covariance of [thetaHat;log(sigmaHat)] and 
% [etaHat,log(sigmaHat)]. 
    methods (Access=protected)  
        
        function C = covThetaHatLogSigmaHat(slme)
%covThetaHatLogSigmaHat - Estimated covariance of [thetaHat;log(sigmaHat)]. 
%   C = covThetaHatLogSigmaHat(slme) returns the estimated covariance
%   matrix of [thetaHat;log(sigmaHat)].

            % (1) Define x = [slme.thetaHat;log(slme.sigmaHat)].
            x = [slme.thetaHat;log(slme.sigmaHat)];

            % (2) Get beta profiled/restricted log likelihood as required.
            switch lower(slme.FitMethod)                                            
               case 'ml'
                   % (2) Get Hessian of negative profiled log likelihood as
                   % a function of [theta;log(sigma)] with beta profiled
                   % out at x = [slme.thetaHat;log(slme.sigmaHat)].
                   fun = makeNegBetaProfiledLogLikelihoodAsFunctionOfThetaLogSigma(slme);                  
               case 'reml'               
                   % (2) Get Hessian of negative restricted log likelihood
                   % as a function of [theta;log(sigma)] at x =
                   % [slme.thetaHat;log(slme.sigmaHat)].                   
                   fun = makeNegRestrictedLogLikelihoodAsFunctionOfThetaLogSigma(slme);                   
            end
            
            % (3) Get Hessian of fun at x. 
            wantRegularized = false;
            H = slme.getHessian(fun,x,wantRegularized);  
            
            % (4) Invert H to get the covariance matrix. H may be singular. 
            % Turn warning MATLAB:singularMatrix off, do the inversion and 
            % restore the warning state.
            warnState = warning('query','all');
            warning('off','MATLAB:singularMatrix');
            warning('off','MATLAB:illConditionedMatrix');
            warning('off','MATLAB:nearlySingularMatrix');
            cleanupObj = onCleanup(@() warning(warnState));
            
            % (5) Invert H. If sigma is fixed, account for that by
            % reporting a 0 variance in the last row and last column of C.
            if slme.isSigmaFixed
                n = size(H,1);
                if n == 1
                    % sigma is the only covariance parameter and it is
                    % fixed. So set it's variance to 0.
                    C = 0;
                else
                    % The last row and last column of H should be 0. Invert
                    % H(1:end-1,1:end-1) and store in C(1:end-1,1:end-1).
                    % The last row and last column of C should be 0 since
                    % sigma is fixed.
                    C = zeros(n);
                    m = n-1;
                    C(1:m,1:m) = H(1:m,1:m) \ eye(m);
                end
            else
                C = H \ eye(size(H));               
                %C = pinv(H);            
            end
                       
        end % covThetaHatLogSigmaHat.        
                
        function C = covEtaHatLogSigmaHat(slme)
%covEtaHatLogSigmaHat - Estimated covariance of [EtaHat;log(sigmaHat)]. 
%   C = covEtaHatLogSigmaHat(slme) returns the estimated covariance
%   matrix of [etaHat;log(sigmaHat)] where etaHat is the estimated Natural
%   parameter vector.            
            
            % (1) Get etaHat and sigmaHat.
            etaHat = getNaturalParameters(slme.Psi);
            sigmaHat = slme.sigmaHat;

            % (2) Define x = [etaHat;log(sigmaHat)].
            x = [etaHat;log(sigmaHat)];

            % (3) Get beta profiled/restricted log likelihood as required.
            switch lower(slme.FitMethod)                                            
               case 'ml'
                   % (3) Get Hessian of negative profiled log likelihood as
                   % a function of [eta;log(sigma)] with beta profiled
                   % out at x.
                   fun = makeNegBetaProfiledLogLikelihoodAsFunctionOfEtaLogSigma(slme);                  
               case 'reml'               
                   % (3) Get Hessian of negative restricted log likelihood
                   % as a function of [eta;log(sigma)] at x.                   
                   fun = makeNegRestrictedLogLikelihoodAsFunctionOfEtaLogSigma(slme);                   
            end
            
            % (4) Get Hessian of fun at x.
            wantRegularized = false;
            try
                H = slme.getHessian(fun,x,wantRegularized);
            catch ME %#ok<NASGU>
                H = NaN(length(x));
            end
            
            % (5) Invert H to get the covariance matrix. H may be singular. 
            % Turn warning MATLAB:singularMatrix off, do the inversion and 
            % restore the warning state.
            warnState = warning('query','all');
            warning('off','MATLAB:singularMatrix');
            warning('off','MATLAB:illConditionedMatrix');
            warning('off','MATLAB:nearlySingularMatrix');
            cleanupObj = onCleanup(@() warning(warnState));
            
            try
                C = slme.covarianceOnNaturalScale(H);
                %C = covarianceOnNaturalScale(slme,H);
                %C = H \ eye(size(H));
            catch ME %#ok<NASGU>
                C = H \ eye(size(H));
            end
            %C = pinv(H);
            
            % (6) If residual std is fixed, the last row and last col of C
            % should be 0.
            if slme.isSigmaFixed
                C(end,:) = 0;
                C(:,end) = 0;
            end
            
        end % end of covEtaHatLogSigmaHat.   
        
%         function C = covarianceOnNaturalScale(slme,H) 
% %covarianceOnNaturalScale - Approximate the covariance of natural parameters.
% %   C = covarianceOnNaturalScale(H) takes a symmetric matrix H where H is
% %   either the Hessian of the negative profiled log likelihood (beta
% %   profiled out) or the Hessian of the negative restricted log likelihood
% %   as a function of [eta;log(sigma)] evaluated at [etaHat;log(sigmaHat)].
% %   The covariance C is given by inv(H) except that here we account for
% %   cases when H is not positive definite.
% 
%             % (1) Tentatively, set C as the pseudo inverse of H.
%             C = pinv(H);
% 
%             % (2) Which elements of the Natural parameter vector are
%             % estimable considering H as equivalent to the Normal matrix
%             % X'*X in linear regression.
%             sizeH = size(H,1);
%             F = internal.stats.isEstimable(eye(sizeH),'NormalMatrix',H,'TolSVD',sqrt(eps(class(H))));
% 
%             % (3) If C(i,i) is < 0, F(i) must be set to false since the ith 
%             % element of Natural parameter vector cannot have a negative
%             % variance.
%             diagC = diag(C);
%             F(diagC < 0) = false;
%             
%             % (4) Set elements of non-estimable rows and columns as NaN.
%             C(~F,:) = NaN;
%             C(:,~F) = NaN;                        
%             
%         end % end of covarianceOnNaturalScale.        

        
    end
              
% Public methods related to fitting.
    methods (Access=public)            
        
        function slme = StandardLinearMixedModel(X,y,Z,Psi,FitMethod,dofit,dostats,varargin)
%

%StandardLinearMixedModel - Create a LME model in standard form.
%   slme = StandardLinearMixedModel(X,y,Z,Psi,FitMethod,dofit,dostats)
%   takes the N by p fixed effects design matrix X, N by 1 response vector
%   y, N by q random effects design matrix Z, an object Psi of type
%   CovarianceMatrix encapsulating the random effects covariance matrix
%   (sigma^2 * D), a string FitMethod (either 'ML' or 'REML') and returns a
%   fitted StandardLinearMixedModel object slme if dofit is true. The
%   residual error is assumed to have an isotropic diagonal covariance with
%   variance sigma^2. If dofit is false then an unfitted object is returned
%   which you must then fit by calling the refit method. If dostats is true
%   then the returned object is ready for doing stats otherwise you must
%   call the initstats method on the object to make it ready for doing
%   stats. If dofit is false then an unfitted object is returned and
%   dostats has no effect.
%
%   In most cases, you would set both dofit and dostats to true. If you
%   just want to access the estimated parameters thetaHat, betaHat, bHat,
%   sigmaHat, loglikHat, Psi then set dofit to true and dostats to false.
%   The idea is that computing stats may not be needed for certain
%   operations such as likelihood ratio tests.
%
%   Before calling any of the stats methods, both isFitToData and
%   isReadyForStats flags must be true.
%
%   slme = StandardLinearMixedModel(...,Name,Value,...) also supplies
%   optional name/value pairs:
%
%           Name                      Value
%           'Optimizer'               A string specifying the name of the 
%                                     optimizer. Supported values are
%                                     'quasinewton' and 'fminunc'.
%
%           'OptimizerOptions'        A structure containing the 
%                                     optimization options to be passed to 
%                                     the optimizer.
%
%           'InitializationMethod'    A string indicating how the initial
%                                     value of parameters should be chosen 
%                                     to initialize the optimizer.
%
%           'CheckHessian'            Either true or false. If true, we
%                                     perform positive definiteness checks
%                                     on (a) Hessian of the objective
%                                     function at convergence (b) Covariance 
%                                     of [thetaHat;log(sigmaHat)] and (c)
%                                     Covariance of [etaHat;log(sigmaHat)].
%                                     Default is false.            
%
%           'ResidualStd'             A positive scalar indicating the
%                                     fixed value of residual standard 
%                                     deviation. Default is NaN indicating
%                                     that residual standard deviation is
%                                     not fixed but estimated from the
%                                     data.

            % (0) No arg constructor
            if ( nargin == 0 )
                return;
            end

            % (1) Ensure that dofit and dostats are scalar logicals.
            assert( isscalar(dofit) & islogical(dofit) );
            assert( isscalar(dostats) & islogical(dostats) );
            
            % (2) Validate inputs.
            %[X,y,Z,Psi,FitMethod] = validateInputs(slme,X,y,Z,Psi,FitMethod);
            
            % (3) Set N, p and q.
            [N,p] = size(X); %#ok<*PROP>
            q = size(Z,2);            
            slme.N = N;
            slme.p = p;
            slme.q = q;
            
            % (4) Set X, y, Z, Psi and FitMethod.
            slme.X = X;
            slme.y = y;
            slme.Z = Z;
            slme.Psi = Psi;
            slme.FitMethod = FitMethod;            
            
            % (5) Parse optimization options.
                % (5a) Default optimization options.
                dfltOptimizer = 'quasinewton';
                dfltOptimizerOptions = struct([]);
                dfltInitializationMethod = 'default';
                dfltCheckHessian = false;
                dfltResidualStd = NaN;
                
                % (5b) Optional parameter names and their default values.
                paramNames = {  'Optimizer',   'OptimizerOptions',   'InitializationMethod',   'CheckHessian',   'ResidualStd'};
                paramDflts = {dfltOptimizer, dfltOptimizerOptions, dfltInitializationMethod, dfltCheckHessian, dfltResidualStd};
           
                % (5c) Parse optional parameter name/value pairs.
                [optimizer,optimizeroptions,initializationmethod,checkhessian,sigmafixed] ...
                    = internal.stats.parseArgs(paramNames,paramDflts,varargin{:});
            
                % 5(d) Validate optimization options.
                [optimizer,optimizeroptions,initializationmethod] ...
                    = validateOptimizationOptions(slme,optimizer,optimizeroptions,initializationmethod);
                
                % 5(e) Validate checkhessian and sigmafixed.
                checkhessian = internal.stats.parseOnOff(checkhessian,'CheckHessian');
                if ~isnan(sigmafixed)
                    sigmafixed = validateSigma(slme,sigmafixed);                    
                end
                    
            % (6) Set optimization options in the object. Also store the
            % Hessian check option selected by the user.
            slme.Optimizer = optimizer;
            slme.OptimizerOptions = optimizeroptions;
            slme.InitializationMethod = initializationmethod;
            slme.CheckHessian = checkhessian;   
            if ~isnan(sigmafixed)
                slme.isSigmaFixed = true;
                slme.sigmaFixed = sigmafixed;
            end
            
            % (7) Fit the standard LME model if asked.
            if (dofit == true)
                slme = refit(slme);
                % (8) Make the object ready for doing stats if required.
                if (dostats == true)
                    slme = initstats(slme);
                end
            end                        
            
        end % end of StandardLinearMixedModel.        
        
        function slme = refit(slme)
%refit - Fit and refit a Linear Mixed Effects (LME) model in standard form.
%   slme = refit(slme) fits a standard LME model previously initialized
%   with the StandardLinearMixedModel method. The public properties of slme
%   such as y, X, etc. can be changed if desired and the model can be refit
%   using the refit method. When changing X, y, Z, Psi or FitMethod, a call
%   to refit and initstats is necessary before calling any of the methods
%   that compute stats on the fitted model.                        

            % (0) Precompute the various matrix and vector products for the
            % current X, y and Z for later use in solving mixed model
            % equations.
            X = slme.X;
            y = slme.y;
            Z = slme.Z;
            slme.XtX = X'*X;
            slme.Xty = X'*y;
            slme.XtZ = X'*Z;
            slme.Zty = Z'*y;
            slme.ZtZ = Z'*Z;

            % (1) Fit LME to get optimal theta and store it in the object.
            slme.thetaHat = solveForThetaHat(slme); 
            
            % (2) Find betaHat and bHat by solving the mixed model eqns. at 
            % slme.thetaHat. Then use betaHat, bHat and slme.thetaHat to 
            % compute sigmaHat and store the results in slme.            
            [slme.betaHat,slme.bHat,slme.sigmaHat,~,~,~,slme.DeltabHat] = ...
                solveMixedModelEquations(slme,slme.thetaHat);            
            
            % (3) Set maximized log likelihood / restricted log likelihood.
            switch lower(slme.FitMethod)
                case 'ml'
                    slme.loglikHat = ...
                        BetaSigmaProfiledLogLikelihood(slme,slme.thetaHat);
                case 'reml'
                    slme.loglikHat = ...
                        SigmaProfiledRestrictedLogLikelihood(slme,slme.thetaHat);
                otherwise
                    % <entry key="BadFitMethod">''FitMethod'' must be either ''ML'' or ''REML''.</entry>
                    error(message('stats:classreg:regr:lmeutils:StandardLinearMixedModel:BadFitMethod'));
            end
            
            % (4) Update slme.Psi to contain optimized theta and sigma.
            slme.Psi = setUnconstrainedParameters(slme.Psi,slme.thetaHat);
            slme.Psi = setSigma(slme.Psi,slme.sigmaHat);                                                
            
            % (5) Mark the object as Fitted.
            slme.isFitToData = true;
            
        end % end of refit.        
        
        function slme = initstats(slme)
%initstats - Make a fitted StandardLinearMixedModel object ready for stats.
%   slme = initstats(slme) takes a fitted StandardLinearMixedModel object
%   slme and makes it ready for stats. slme must be a fitted object
%   returned by either StandardLinearMixedModel constructor or the refit
%   method.
            
            % (1) Ensure that slme has been fit to data.
            if (slme.isFitToData == false)
                % <entry key="MustRefitFirst">Call the method refit first and then initstats.</entry>
                error(message('stats:classreg:regr:lmeutils:StandardLinearMixedModel:MustRefitFirst'));
            end

            % (2) Check the Hessian of objective function at convergence.
            if slme.CheckHessian == true
                checkObjectiveFunctionHessianAtThetaHat(slme);
            end
            
            % (3) Store the rank of X. Since we have already ensured that X
            % is full column rank, don't compute the rank again.
            % slme.rankX = rank(slme.X);
            slme.rankX = slme.p;

            % (4) Store estimated covariance of betaHat in slme.
            slme.covbetaHat = ...
                covBetaHatAsFunctionOfThetaSigma(slme,slme.thetaHat,slme.sigmaHat);
            
            % (5) Store estimated covariance of [thetaHat;log(sigmaHat)] in
            % slme.
            slme.covthetaHatlogsigmaHat = covThetaHatLogSigmaHat(slme);

            % (6) Ensure that slme.covthetaHatlogsigmaHat is positive
            % definite.
            if slme.CheckHessian == true
                % <entry key="Message_NotSPDCovarianceUnconstrainedScale">The covariance matrix of covariance parameters on the unconstrained scale is not positive definite. This may indicate a model with more covariance parameters than those supported by data.</entry>
                msg1ID = 'stats:classreg:regr:lmeutils:StandardLinearMixedModel:Message_NotSPDCovarianceUnconstrainedScale';
                slme.checkPositiveDefinite(slme.covthetaHatlogsigmaHat,msg1ID);
            end
            
            % (7) Store estimated covariance of [etaHat;log(sigmaHat)] in
            % slme.
            slme.covetaHatlogsigmaHat = covEtaHatLogSigmaHat(slme);
            
            % (8) Ensure that slme.covetaHatlogsigmaHat is positive
            % definite.
            if slme.CheckHessian == true
                % <entry key="Message_NotSPDCovarianceNaturalScale">The covariance matrix of covariance parameters on the Natural scale is not positive definite. This may indicate a model with more covariance parameters than those supported by data.</entry>
                msg2ID = 'stats:classreg:regr:lmeutils:StandardLinearMixedModel:Message_NotSPDCovarianceNaturalScale';
                slme.checkPositiveDefinite(slme.covetaHatlogsigmaHat,msg2ID);
            end
            
            % TODO: This matrix may be very large. Decide whether to save
            % it in the object or compute it on the fly.
            % (9) Store estimated covariance of [betaHat-beta;bHat-b] in
            % slme.
            %slme.covbetaHatbHat = ...
            %    covBetaHatBHatAsFunctionOfThetaSigma(slme,slme.thetaHat,slme.sigmaHat);                        
            
            % (10) Mark the object as ready for stats.
            slme.isReadyForStats = true;            
            
        end % end of initstats.                
        
    end           
    
% Public methods to get Satterthwaite DF for T/F tests on beta/b vector.
    methods (Access=public)      
        
        function df = dfBetaTTest(slme,c)
%dfBetaTTest - Satterthwaite degrees of freedom for testing c'*beta = e.
%   df = dfBetaTTest(slme,c) returns the Satterthwaite denominator degrees
%   of freedom for the test: c'*beta = e. beta is a p by 1 vector and so c 
%   must be a p by 1 vector as well. The calculation is insensitive to the 
%   value of e.
    
            % (1) Ensure that c has the right size.
            if size(c,1) == 1
                c = c';
            end
            assert( all(size(c) == [slme.p,1]) );
            
            % (2) If V denotes the covariance of betaHat at the values 
            % (theta,sigma) then make the function c'*V(theta,sigma)*c as a
            % function of [theta;log(sigma)].
            ctVcfun = makecovcTBetaHatAsFunctionOfThetaLogSigma(slme,c);            
            
            % (3) Let x = [slme.thetaHat;log(slme.sigmaHat)].
            x = [slme.thetaHat;log(slme.sigmaHat)];
            
            % (4) gHat is the gradient of ctVcfun evaluated at x.
            gHat = slme.getGradient(ctVcfun,x);
            
            % (5) Get CHat the covariance matrix of [thetaHat;log(sigmaHat)]
            % at x. This has been already computed by refit.
            CHat = slme.covthetaHatlogsigmaHat;
            
            % (6) Compute df.
            df = 2*(ctVcfun(x))^2/(gHat'*CHat*gHat);            
            
            % (7) Ensure that df is never < 0.
            df = max(0,df);
            
        end % end dfBetaTTest.                        
        
        function df = dfBetaFTest(slme,L)
%dfBetaFTest - Satterthwaite degrees of freedom for testing L*beta = e.           
%   df = dfBetaFTest(slme,L) returns the Satterthwaite denominator degrees 
%   of freedom for the test L*beta = e. beta is the p by 1 vector and L 
%   must be a r by p matrix with r <= p and rank(L) = r. The calculation is
%   insensitive to the value of e.
            
           % (1) Make sure L is sensible.
           r = size(L,1);
           assert( size(L,2) == slme.p & r <= slme.p );
           assert( rank(L) == r );

           % (2) Get the estimated covariance of betaHat.
           V = slme.covbetaHat;
           
           % (3) Get eigendecomposition of L*V*L' = Veig*Lambdaeig*Veig'
           [Veig,Lambdaeig] = eig(L*V*L');
           Lambdaeig = diag(Lambdaeig);
           
           % (4) Get U matrix.
           U = Veig';
           
           % (5) Get B such that B(i,:) = b_i^T.
           B = U*L;
           
           % (6) Get [thetaHat;log(sigmaHat)].
           x = [slme.thetaHat;log(slme.sigmaHat)];
           
           % (7) Get covariance of [thetaHat;log(sigmaHat)].
           C = slme.covthetaHatlogsigmaHat;
           
           % (8) Initialize nu and loop.
           nu = zeros(r,1);           
           for i = 1:r
                % (9) Get covariance of b_i'*betaHat as a function of
                % [theta;log(sigma)].
                bi = B(i,:)';
                fi = makecovcTBetaHatAsFunctionOfThetaLogSigma(slme,bi);
                
                % (10) Get gradient of fi at x = [thetaHat;log(sigmaHat)].
                gi = slme.getGradient(fi,x);
                
                % (11) Get nu(i)
                nu(i) = 2*(Lambdaeig(i)^2)/(gi'*C*gi);                
           end
           
           % (12) Compute G.
           nu = nu( (nu > 2) );
           if isempty(nu)
                G = 0;
           else
                G = sum( (nu./(nu-2)) );
           end
           %G = sum( (nu./(nu-2)) .* (nu > 2) );
           
           % (13) Compute df.
           if ( G > r )
                df = 2*G/(G-r);            
           else
                df = 0;
           end
           
        end % end of dfBetaFTest.
        
        function df = dfBetaBTTest(slme,c)
%dfBetaBTTest - Satterthwaite degrees of freedom for testing c'*[beta;b] = e.
%   df = dfBetaBTTest(slme,c) returns the Satterthwaite denominator degrees
%   of freedom for the test: c'*[beta;b] = e. beta is a p by 1 vector and b
%   is a q by 1 vector and c must be a (p+q) by 1 vector. The calculation
%   is insensitive to the value of e.
            
            % (1) Ensure that c has the right size.
            if size(c,1) == 1
                c = c';
            end
            assert( all(size(c) == [slme.p+slme.q,1]) );
            
            % (2) If V denotes the covariance of [betaHat-beta;bHat-b] at
            % the values (theta,sigma) then make the function
            % c'*V(theta,sigma)*c as a function of [theta;log(sigma)].
            ctVcfun = makecovcTBetaHatBHatAsFunctionOfThetaLogSigma(slme,c);
            
            % (3) Let x = [slme.thetaHat;log(slme.sigmaHat)].
            x = [slme.thetaHat;log(slme.sigmaHat)];
            
            % (4) gHat is the gradient of ctVcfun evaluated at x.
            gHat = slme.getGradient(ctVcfun,x);
            
            % (5) Get CHat the covariance matrix of [thetaHat;log(sigmaHat)]
            % at x. This has been already computed by refit.
            CHat = slme.covthetaHatlogsigmaHat;
            
            % (6) Compute df.
            df = 2*(ctVcfun(x))^2/(gHat'*CHat*gHat);
            
            % (7) Ensure that df is never < 0.
            df = max(0,df);
            
        end %end of dfBetaBTTest.

        function df = dfBetaBFTest(slme,L)
%dfBetaBFTest - Satterthwaite degrees of freedom for testing L*[beta;b] = e.           
%   df = dfBetaBFTest(slme,L) returns the Satterthwaite denominator degrees 
%   of freedom for the test L*[beta;b] = e. beta is the p by 1 vector, b is
%   a q by 1 vector and L is a r by (p+q) matrix with r <= (p+q) and with 
%   rank(L) = r. The calculation is insensitive to the value of e.            
            
           % (1) Make sure L is sensible.
           r = size(L,1);
           assert( size(L,2) == (slme.p+slme.q) & r <= (slme.p+slme.q) );
           assert( rank(L) == r );
                      
           % (2) Get L*V*L'. Here's one way to do this:
           %    (a) Get the estimated covariance of [betaHat-beta;bHat - b].
           %        V = covBetaHatBHatAsFunctionOfThetaSigma(slme,slme.thetaHat,slme.sigmaHat);
           %    (b) Form L*V*L'
           % What follows is a more efficient way of doing the same thing.
           p = slme.p;
           q = slme.q;
           Xnew = L(:,1:p);
           Znew = L(:,p+1:p+q);
           [~,LVLt] = getEstimateAndVariance(slme,Xnew,Znew,slme.betaHat,...
               slme.DeltabHat,slme.thetaHat,slme.sigmaHat,'covariance');
           
           % (3) Get eigendecomposition of L*V*L' = Veig*Lambdaeig*Veig'
           % [Veig,Lambdaeig] = eig(L*V*L');
           [Veig,Lambdaeig] = eig(LVLt);
           Lambdaeig = diag(Lambdaeig);
           
           % (4) Get U matrix.
           U = Veig';
           
           % (5) Get B such that B(i,:) = b_i^T.
           B = U*L;
           
           % (6) Get [thetaHat;log(sigmaHat)].
           x = [slme.thetaHat;log(slme.sigmaHat)];
           
           % (7) Get covariance of [thetaHat;log(sigmaHat)].
           C = slme.covthetaHatlogsigmaHat;
           
           % (8) Initialize nu and loop.
           nu = zeros(r,1);           
           for i = 1:r
                % (9) Get covariance of b_i'*betaHat as a function of
                % [theta;log(sigma)].
                bi = B(i,:)';
                fi = makecovcTBetaHatBHatAsFunctionOfThetaLogSigma(slme,bi);
                % (10) Get gradient of fi at x = [thetaHat;log(sigmaHat)].
                gi = slme.getGradient(fi,x);
                
                % (11) Get nu(i)
                nu(i) = 2*(Lambdaeig(i)^2)/(gi'*C*gi);                
           end
           
           % (12) Compute G.
           G = sum( (nu./(nu-2)) .* (nu > 2) );
           
           % (13) Compute df.
           if ( G > r )
                df = 2*G/(G-r);            
           else
                df = 0;
           end

        end % end dfBetaBFTest.
        
    end
                     
% Public method for computing model criterion.
    methods (Access=public)     
        
        function crittable= modelCriterion(slme)
%modelCriterion - Compute table containing model criterion info.
%   crittable= modelCriterion(slme) takes a StandardLinearMixedModel object
%   slme and computes a table of model criterion such as AIC, BIC, logLik 
%   and Deviance.

    % Create a structure stats for classreg.regr.modelutils.modelcriterion
    % such that:
    % L = stats.LogLikelihood;
    % k = stats.NumCoefficients;
    % n = stats.NumObservations;      
    
            % (1) Get N and p.
            N = slme.N;
            p = slme.p;
            
            % (2) Total number of unknown parameters in the model.
            if slme.isSigmaFixed
                % residual std. was not estimated.
                stats.NumCoefficients = ...
                slme.Psi.NumParametersExcludingSigma + p;
            else
                % residual std. was estimated.
                stats.NumCoefficients = ...
                    slme.Psi.NumParametersExcludingSigma + (p+1);
            end
            
            % (3) Get effective number of observations based on FitMethod.
            switch lower(slme.FitMethod)
                case {'ml'}
                    stats.NumObservations = N;
                case {'reml'}
                    stats.NumObservations = (N-p);
                otherwise
                    error(message('stats:classreg:regr:lmeutils:StandardLinearMixedModel:BadFitMethod'));                    
            end
            
            % (4) Set maximized log likelihood.
            loglikHat = slme.loglikHat;
            stats.LogLikelihood = loglikHat;
            
            % (5) Call modelcriterion utility function.
            crit = classreg.regr.modelutils.modelcriterion(stats,'all',true);
            
            % (6) Get deviance and create output table.
            Deviance = -2*loglikHat;
            crittable = table(crit.AIC, crit.BIC, loglikHat, Deviance,...
                'VariableNames',{'AIC' 'BIC' 'logLik' 'Deviance'});
            
        end % end of modelCriterion.        
        
    end
           
% Public method to get estimated posterior covariance of random effects vector.
    methods (Access=public)      
        
        function C = postCovb(slme)
%postCovb - Posterior covariance of b. 
%   C = postCovb(slme) takes a StandardLinearMixedModel object slme and
%   computes the estimated posterior covariance matrix of b.

            % The posterior covariance is given by:
            % C = sigma^2 * inv(Z'*Z + Delta'*Delta) 
            % with Delta'*Delta = inv(D).
            % We replace sigma by sigmaHat.

            % (1) Get Lambda = inv(Delta).
            Lambda = getLowerTriangularCholeskyFactor(slme.Psi);
            
            % (2) Get Z and q.
            Z = slme.Z;
            q = slme.q;
            
            % (3) Form the matrix U.
            U = Z*Lambda;
            
            % (4) Get R and S such that S*R'*R*S' = U'*U + eye(q) where R
            % and S are q by q. S is a permutation matrix: S*S' = eye(q).
            Iq = spdiags(ones(q,1),0,q,q);
            [R,status,S] = chol(U'*U + Iq);
            
            % (5) Ensure that Cholesky factorization worked.
            if ( status ~= 0 )
                % <entry key="ErrorPostCovb">An error occured in getting the posterior covariance of random effects vector.</entry>
                error(message('stats:classreg:regr:lmeutils:StandardLinearMixedModel:ErrorPostCovb'));
            end
            
            % (6) Get sigmaHat.
            sigmaHat = slme.sigmaHat;
            
            % (7) C can also be written as:
            %     C = sigma^2 * (Lambda*S*inv(R) * inv(R)'*S'*Lambda').
            % We replace sigma by sigmaHat.
            T = (Lambda*S) / R;
            C = (sigmaHat^2) * (T*T');
            
        end % end of postCovb.    
        
    end
    
% Public methods related to fitted values.
    methods (Access=public) 
        
        function yfit = fitted(slme,wantConditional)
%fitted - Returns the fitted response from a StandardLinearMixedModel.
%   yfit = fitted(slme,true) returns the N by 1 vector representing the
%   fitted conditional response from the StandardLinearMixedModel slme such
%   that yfit = X*betaHat + Z*bHat.
%
%   yfit = fitted(slme,false) returns the N by 1 vector representing the
%   fitted marginal response from the StandardLinearMixedModel slme such
%   that yfit = X*betaHat.

            % (1) Ensure that wantConditional is sensible.
            assert( islogical(wantConditional) & isscalar(wantConditional) );

            % (2) Decide which type of fitted values to return.
            if ( wantConditional == true )
                yfit = slme.X*slme.betaHat + slme.Z*slme.bHat;
            else
                yfit = slme.X*slme.betaHat;
            end

        end % end of fitted.       
        
    end
 
% Public method to make predictions on new data.
    methods (Access=public)   
        
        function [ypred,CI,DF] = predict(slme,Xnew,Znew,alpha,dfmethod,wantConditional,wantPointwise,wantCurve,hasIntercept)
%predict - Predictions from a fitted StandardLinearMixedModel.
%   [ypred,CI,DF] = predict(slme,Xnew,Znew,alpha,dfmethod,wantConditional,wantPointwise,wantCurve,hasIntercept)
%   takes a fitted StandardLinearMixedModel slme, a fixed effects design
%   matrix Xnew, a random effects design matrix Znew, a confidence level
%   alpha, a degrees of freedom method dfmethod, flags for conditional,
%   pointwise and curve predictions and outputs predictions and (1-alpha)
%   confidence intervals (CIs) for the true predictions from slme at Xnew
%   and Znew. If wantConditional is true then the predictions include
%   contributions from both fixed effects and BLUPs of random effects
%   otherwise the predictions include contributions from only the fixed
%   effects. If wantPointwise is true then pointwise CIs are computed
%   otherwise simultaneous CIs are computed. If wantCurve is true then CIs
%   are for the curve excluding variability due to observation noise
%   otherwise variability due to observation noise is included. Xnew must
%   be M-by-p and Znew must be M-by-q where p = slme.p and q = slme.q. The
%   output ypred is a M-by-1 vector containing the predictions
%   corresponding to the rows of Xnew and Znew and output CI is a M-by-2
%   matrix. Each row of CI is a (1-alpha) confidence interval for the true
%   prediction such that CI(:,1) contains the lower confidence limit and
%   CI(:,2) contains the upper confidence limit. DF is a M-by-1 vector
%   containing the DF values used in computing the CIs if wantPointwise is
%   true. If wantPointwise is false then DF is a scalar containing the DF
%   value used for computing the Scheffe simultaneous CIs. hasIntercept
%   should be set to true if the fixed effects part of the model has an
%   intercept and false otherwise.

            args = {Xnew,Znew,alpha,dfmethod,wantConditional,wantPointwise,wantCurve,hasIntercept};
            switch nargout
                case {0,1}
                    ypred         = predict@classreg.regr.lmeutils.StandardLinearLikeMixedModel(slme,args{:});
                case 2
                    [ypred,CI]    = predict@classreg.regr.lmeutils.StandardLinearLikeMixedModel(slme,args{:});
                case 3
                    [ypred,CI,DF] = predict@classreg.regr.lmeutils.StandardLinearLikeMixedModel(slme,args{:});
            end
            %[ypred,CI,DF] = predict@classreg.regr.lmeutils.StandardLinearLikeMixedModel(slme,...
            %    Xnew,Znew,alpha,dfmethod,wantConditional,wantPointwise,wantCurve,hasIntercept);

        end % end of predict.
        
    end
    
% Public method for generating random data from the fitted model.
    methods (Access=public)  
        
        function ysim = random(slme,S,Xsim,Zsim)
%random - Generate random data from a fitted StandardLinearMixedModel.
%   ysim = random(slme,S,Xsim,Zsim) takes a fitted StandardLinearMixedModel
%   slme and generates random data from slme using the RandStream object S
%   at the specified fixed effects design matrix Xsim and random effects
%   design matrix Zsim. Xsim must be M by p and Zsim must be M by q where p
%   = slme.p and q = slme.q. The output ysim is a M by 1 vector. Set S = []
%   to use the default RandStream object.
    
            % (1) Ensure that Xsim is sensible.
            assert( isnumeric(Xsim) & isreal(Xsim) & ismatrix(Xsim) );
            assert( size(Xsim,2) == slme.p );
            
            % (2) Ensure that Zsim is sensible.
            assert( isnumeric(Zsim) & isreal(Zsim) & ismatrix(Zsim) );
            assert( size(Zsim,2) == slme.q );
            
            % (3) Ensure that Xsim and Zsim have the same number of rows.
            assert( size(Xsim,1) == size(Zsim,1) );
            
            % (4) Extract Psi, betaHat, sigmaHat from slme. 
            Psi = slme.Psi;
            betaHat = slme.betaHat;
            sigmaHat = slme.sigmaHat;
            
            % (5) Also extract NumReps from Psi.
            NumReps = Psi.NumReps;

            % (6) Psi is a BlockedCovariance object. It has NumBlocks smaller
            % covariance matrices given in the Psi.Matrices cell array.
            % Psi.Matrices{i} is repeated NumReps(i) times along the block
            % diagonal of Psi. Psi.Matrices{i} has size SizeVec(i).
            
            % (7) The combined random effects vector bsim will be of size 
            % sum_{i=1 to NumBlocks} SizeVec(i)*NumReps(i).
            bsim = randomb(slme,S,NumReps);            
                       
            % (8) Generate the noise term epsilon ~ N(0,sigma^2 I_N). 
            N = size(Xsim,1);
            if isempty(S)
                % Use default RandStream.
                epsilonsim = sigmaHat*randn(N,1);
            else
                % Use user supplied RandStream S.
                epsilonsim = sigmaHat*randn(S,N,1);
            end
            
            % (9) Generate simulated data.
            if isempty(bsim)
                % Only fixed effects.
                ysim = Xsim * betaHat + epsilonsim;
            else
                ysim = Xsim * betaHat + Zsim * bsim + epsilonsim;
            end

        end % end of random.        
        
    end
    
% Protected methods to get raw/Pearson/Standardized residuals.    
    methods (Access=protected)      
        
        % Raw residuals (marginal and conditional).
        function rawr = getRawResiduals(slme,wantConditional)
%getRawResiduals - Get Raw residuals.
%   rawr = getRawResiduals(slme,wantConditional) takes an object slme of
%   type StandardLinearMixedModel and returns the raw conditional residuals
%   if wantConditional is true and the raw marginal residuals if the flag
%   wantConditional is false.
           
            % (1) Ensure that wantConditional is sensible.
            assert( isscalar(wantConditional) & islogical(wantConditional) );
            
            % (2) Decide which type of raw residuals to return.
            if ( wantConditional == true )
                % Conditional residuals. 
                rawr = slme.y - slme.X*slme.betaHat - slme.Z*slme.bHat; 
            else
                % Marginal residuals. 
                rawr = slme.y - slme.X*slme.betaHat;
            end

        end % end of getRawResiduals.
        
        % Pearson residuals (marginal and conditional).
        function pearsonr = getPearsonResiduals(slme,wantConditional)
%getPearsonResiduals - Get Pearson residuals.           
%   pearsonr = getPearsonResiduals(slme,wantConditional) takes an object
%   slme of type StandardLinearMixedModel and returns the Pearson
%   conditional residuals if wantConditional is true and the Pearson
%   marginal residuals if wantConditional is false.

            % (1) Ensure that wantConditional is sensible.
            assert( isscalar(wantConditional) & islogical(wantConditional) );

            % (2) Decide which type of Pearson residuals to return.
            if ( wantConditional == true )
                % Conditional residuals.
                
                % (1) Get raw conditional residuals.
                rawr = getRawResiduals(slme,true);
                
                % (2) Divide rawr elementwise by a vector with all elements
                % equal to slme.sigmaHat: stdvec = slme.sigmaHat*ones(slme.N,1);
                stdvec = slme.sigmaHat;
                
                % (3) Form the conditional, Pearson residuals.
                pearsonr = rawr ./ stdvec;
                
            else
                % Marginal residuals.
                
                % (1) Get raw marginal residuals.
                rawr = getRawResiduals(slme,false);
                
                % Divide rawr elementwise by the sqrt of the diagonal of M
                % where:
                %
                % M = sigma^2 * inv(H) where H is a matrix such that
                %
                % inv(H) = Z*D*Z' + IN and IN = eye(N). 
                %
                % If D = Lambda*Lambda' and U = Z*Lambda then we can write:
                % 
                % inv(H) = U*U' + IN. 
                %
                % So the M matrix is: 
                % 
                % M = sigma^2 * inv(H) = sigma^2 * (U*U' + IN)
                
                % (2) Get sigmaHat and Lambda.
                sigmaHat = slme.sigmaHat;
                Lambda = getLowerTriangularCholeskyFactor(slme.Psi);
                
                % (3) Form the U matrix.
                U = slme.Z * Lambda;    
                
                % (4) Get diagonal of invH and the appropriate std vector. 
                diaginvH = ( sum(U.^2,2) + 1 );
                stdvec = sigmaHat * sqrt(diaginvH);                
                
                % (5) Form the marginal, Pearson residuals.
                pearsonr = rawr ./ stdvec;                
            end
            
        end % end of getPearsonResiduals.        

        % Internally Studentized residuals (marginal and conditional).
        % These are called 'Standardized' residuals in LinearModel.
        function studr = getStudentizedResiduals(slme,wantConditional)
%getStudentizedResiduals - Get Studentized residuals.
%   studr = getStudentizedResiduals(slme,wantConditional) takes an object
%   slme of type StandardLinearMixedModel and returns the Studentized
%   conditional residuals if wantConditional is true and the Studentized
%   marginal residuals if wantConditional is false.

            % (1) Ensure that wantConditional is sensible.
            assert( isscalar(wantConditional) & islogical(wantConditional) );

            % (2) Get sigmaHat, X, Z, q, Lambda and U.
            sigmaHat = slme.sigmaHat;
            X = slme.X;
            Z = slme.Z;
            q = slme.q;
            Lambda = getLowerTriangularCholeskyFactor(slme.Psi);
            U = Z*Lambda;
            
            % (3) Get R and S such that S*R'*R*S' = U'*U + eye(q) where
            % R and S are q by q. S is a permutation matrix: S*S' = eye(q).
            Iq = spdiags(ones(q,1),0,q,q);
            [R,status,S] = chol(U'*U + Iq);
            
            % (4) Ensure that Cholesky factorization worked.
            if ( status ~= 0 )
                % <entry key="ErrorStandardizedResiduals">An error occured in getting standardized residuals.</entry> 
                error(message('stats:classreg:regr:lmeutils:StandardLinearMixedModel:ErrorStandardizedResiduals'));                
            end
            
            % (5) Write inv(Z'*Z + Delta'*Delta) = T*T' where T is such
            % that T = Lambda*S*inv(R).
            T = Lambda * S / R;
            
            % (6) Write
            % X'*H*X = X'*X - X'*Z*inv(Z'*Z + Delta'*Delta)*Z'*X
            %        = X'*X - X'*Z*T*T'*Z'*X
            %        = R1*R1' where R1 is lower triangular.
            try
                Q = X'*Z*T;
                R1 = chol(X'*X - Q*Q','lower');
            catch ME %#ok<NASGU>
                R1 = classreg.regr.lmeutils.covmats.CovarianceMatrix.singularLowerChol(X'*X - Q*Q');
            end
                        
            % (7) Decide which type of Studentized residuals to return.
            if ( wantConditional == true )
                % Conditional residuals.
               
                % (1) Get the raw conditional residuals.
                rawr = getRawResiduals(slme,true);
                
                % We need to divide rawr elementwise by the sqrt of the
                % diagonal of the matrix M where 
                %
                % M = sigma^2 * K * F *K = sigma^2 * H*F*H where
                %
                % H = IN - Z*inv(Z'*Z + Delta'*Delta)*Z' and
                % F = inv(H) - X*inv(X'*H*X)*X'
                %
                % In terms of factors T and R1 (described above), we can
                % write:
                %
                % M = sigma^2 * (H - H*X*inv(X'*H*X)*X'*H)   
                %   = sigma^2 * (IN - Z*inv(Z'*Z + Delta'*Delta)*Z' - H*X*inv(X'*H*X)*X'*H)  
                %   = sigma^2 * (IN - Z*T*T'*Z' - H*X*inv(R1*R1')*X'*H) 
                %   = sigma^2 * (H*F*H)
                                
                % (2) Get H*X.
                HX = X - (Z*T)*(T'*Z'*X);
                
                % (3) Get the diagonal elements of H*F*H.
                diagHFH = 1 - sum((Z*T).^2,2) - sum((HX / R1').^2,2);
                
                % (4) Get the appropriate std vector.
                stdvec = sigmaHat*sqrt(diagHFH);
                
                % (5) Form the conditional Studentized residuals.
                studr = rawr ./ stdvec;                                                                                                
                
            else
                % Marginal residuals.
                
                % (1) Get the raw marginal residuals.
                rawr = getRawResiduals(slme,false);
                
                % We need to divide rawr elementwise by the sqrt of the
                % diagonal of the matrix M where 
                %
                % M = sigma^2 * F   where
                %
                % H = IN - Z*inv(Z'*Z + Delta'*Delta)*Z' and
                % F = inv(H) - X*inv(X'*H*X)*X'
                %
                % inv(H) = Z*D*Z' + IN where D = Lambda*Lambda'
                %
                % In terms of factors U and R1 (described above), we can
                % write:
                %
                % F = U*U' + IN - X*inv(R1*R1')*X'
                
                % (2) Get diagonal of F.
                diagF = 1 + sum((U.^2),2) - sum((X/R1').^2,2);
                
                % (3) Get the appropriate std vector.
                stdvec = sigmaHat*sqrt(diagF);
                
                % (4) Form the marginal Studentized residuals.
                studr = rawr ./ stdvec;                                                
            end
                        
        end % end of getStudentizedResiduals.        
        
    end        
        
end