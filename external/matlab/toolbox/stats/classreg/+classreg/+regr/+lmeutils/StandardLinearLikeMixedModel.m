classdef (Abstract) StandardLinearLikeMixedModel
%  This feature is intended for internal use only and is subject to change
%  at any time without warning. 

%StandardLinearLikeMixedModel - An abstract class to represent aspects that 
%   are common to fitted linear mixed effects (LME) and generalized linear
%   mixed effects (GLME) models in standard form.
%% GLME model.
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
%% LME model.
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

    properties (Abstract=true, GetAccess=public, SetAccess=public)
%y - N-by-1 response vector used to fit the model.
        y

%X - N-by-p fixed effects design matrix used to fit the model.
        X

%Z - N-by-q random effects design matrix used to fit the model.
        Z
        
%Psi - An object of type CovarianceMatrix representing the prior covariance
%   of random effects. Psi = sigma^2 * D(theta). 
        Psi
        
%FitMethod - A string representing the method used to fit the model.
        FitMethod                
    end

    properties (Abstract=true, GetAccess=public, SetAccess=protected)
%N - Number of observations. This is the number of rows in y, X and Z.
        N

%p - Number of fixed effects. This is the number of columns in X.
        p

%q - Number of random effects. This is the number of columns in Z.
        q

%rankX - Rank of the fixed effects design matrix.
        rankX
    end       
    
    properties (Abstract=true, GetAccess=public, SetAccess=protected)
%Optimizer - A string containing the name of the optimizer.
        Optimizer

%OptimizerOptions - A structure containing the optimizer options.
        OptimizerOptions

%CheckHessian - A logical scalar indicating whether to perform the Hessian
%   checks after model fitting.
        CheckHessian
    end

    properties (Abstract=true, GetAccess=public, SetAccess=public)
%InitializationMethod - A string containing the method used to compute the
%   initial values of parameters to start the optimization.
        InitializationMethod
    end
    
    properties (Abstract=true, GetAccess=public, SetAccess=protected)
%betaHat - Estimated fixed effects vector.
        betaHat

%bHat - Estimated random effects vector.
        bHat

%sigmaHat - Estimated residual standard deviation or square root of 
%   dispersion parameter.
        sigmaHat

%thetaHat - Estimated vector of unconstrained parameters for matrix D.
        thetaHat

%loglikHat - Maximized log likelihood or maximized restricted log 
%   likelihood.
        loglikHat
    end

    properties (Abstract=true, GetAccess=public, SetAccess=protected)
%covbetaHat - Estimated covariance of betaHat.
        covbetaHat

%covthetaHatlogsigmaHat - Estimated covariance of [thetaHat;log(sigmaHat)].
        covthetaHatlogsigmaHat

%covetaHatlogsigmaHat - Estimated covariance of [etaHat;log(sigmaHat)].
%   etaHat is the Natural parameter vector corresponding to thetaHat.
        covetaHatlogsigmaHat

%covbetaHatbHat - Estimated covariance of [betaHat-beta;bHat-b] or the
%   joint covariance of [beta;b] given y, thetaHat and sigmaHat.
        covbetaHatbHat
    end
 
    properties (Abstract=true, GetAccess=public, SetAccess=protected)
%isSigmaFixed - True if residual standard deviation sigma is fixed.
        isSigmaFixed
        
%sigmaFixed - Scalar fixed value for residual standard deviation.
        sigmaFixed
    end
    
    properties (Abstract=true, GetAccess=public, SetAccess=protected)
%isFitToData - True if this is a fitted object.
        isFitToData

%isReadyForStats - True if this fitted object is ready for stats.
        isReadyForStats
    end
            
    properties (Constant=true, Hidden=true)               
%AllowedDFMethods - The DF methods that we currently support.        
        AllowedDFMethods = {'None','Residual','Satterthwaite'};  
        
%AllowedResidualTypes - The Residual types that we currently support.
        AllowedResidualTypes = {'Raw','Pearson','Standardized'};
        
%AllowedOptimizers - A list of allowed optimizers.
        AllowedOptimizers = {'fminsearch','fminunc','quasinewton'};     
        
%AllowedInitializationMethods - A list of allowed initialization methods.
        AllowedInitializationMethods = {'random','default'};                
    end
    
    % Public abstract methods related to fitting and stats.    
    methods (Abstract=true, Access=public)   
        slme = refit(slme)
%slme = refit(slme) takes a StandardLinearLikeMixedModel slme and refits 
%   the model. This method should set the following properties: thetaHat, 
%   betaHat, bHat, sigmaHat, loglikHat, Psi and isFitToData.

        slme = initstats(slme)      
%slme = initstats(slme) takes a StandardLinearLikeMixedModel slme and makes 
%   it ready for stats. This method should set the following properties:
%   rankX, covbetaHat, covthetaHatlogsigmaHat, covetaHatlogsigmaHat and 
%   isReadyForStats.
    end
    
    % Public abstract methods related to Satterthwaite DF computation.
    methods (Abstract=true, Access=public)        
        df = dfBetaTTest(slme,c)
%df = dfBetaTTest(slme,c) computes the Satterthwaite degrees of freedom for 
%   testing c'*beta = e.

        df = dfBetaBTTest(slme,c)
%df = dfBetaBTTest(slme,c) computes the Satterthwaite degrees of freedom 
%   for testing c'*[beta;b] = e.

        df = dfBetaFTest(slme,L)
%df = dfBetaFTest(slme,L) computes the Satterthwaite degrees of freedom for 
%   testing L*beta = e.

        df = dfBetaBFTest(slme,L)
%df = dfBetaBFTest(slme,L) computes the Satterthwaite degrees of freedom 
%   for testing L*[beta;b] = e.        
    end        
        
    % Public abstract API methods.
    methods (Abstract=true, Access=public)                        
        crittable= modelCriterion(slme)  
        yfit = fitted(slme,wantConditional)
        ysim = random(slme,S,Xsim,Zsim)
        C = postCovb(slme)
    end
    
    % Protected abstract helper methods.
    methods (Abstract=true, Access=protected)
        r = getRawResiduals(slme,wantConditional)
        r = getPearsonResiduals(slme,wantConditional)
        r = getStudentizedResiduals(slme,wantConditional)        

        C = covEtaHatLogSigmaHat(slme)        
        C = covBetaHatBHat(slme)
        [pred,varpred] = getEstimateAndVariance(slme,Xnew,Znew,betaBar,DeltabBar,theta,sigma)
    end
                                           
    % Protected abstract methods related to validation.
    methods (Abstract=true, Access=protected)
        FitMethod = validateFitMethod(slme,FitMethod)
    end    
    
    % Protected input validation methods.
    methods (Access=protected)       
        function [X,y,Z,Psi,FitMethod] = validateInputs(slme,X,y,Z,Psi,FitMethod)
%validateInputs - Basic checks on inputs. Return the validated inputs.
%   [X,y,Z,Psi,FitMethod] = validateInputs(slme,X,y,Z,Psi,FitMethod) takes
%   user supplied inputs, validates them and returns the validated results.

            % Assert that:
            % (1) X is a numeric, real, matrix of size N by p.
            [N,p] = size(X);
            slme.assertThat(isnumeric(X) & isreal(X) & ismatrix(X),'stats:classreg:regr:lmeutils:StandardLinearLikeMixedModel:BadX');
            
            % (2) X must not have any NaN or Inf values.
            slme.assertThat(~slme.hasNaNInf(X),'stats:classreg:regr:lmeutils:StandardLinearLikeMixedModel:NoNaNInfAllowed','X');
            
            % (3) X has rank p.
            slme.assertThat(rank(X) == p,'stats:classreg:regr:lmeutils:StandardLinearLikeMixedModel:MustBeFullRank_X');                        
            
            % (4) y is a numeric, real, vector of size N by 1.
            slme.assertThat(isnumeric(y) & isreal(y) & isvector(y),'stats:classreg:regr:lmeutils:StandardLinearLikeMixedModel:BadY',num2str(N)); 
            slme.assertThat(       all(size(y) == [N,1])          ,'stats:classreg:regr:lmeutils:StandardLinearLikeMixedModel:BadY',num2str(N)); 

            % (5) y must not have any NaN or Inf values.
            slme.assertThat(~slme.hasNaNInf(y),'stats:classreg:regr:lmeutils:StandardLinearLikeMixedModel:NoNaNInfAllowed','y');
            
            % (6) Z is a numeric, real, matrix of size N by q.
            q = size(Z,2);
            slme.assertThat( isnumeric(Z) & isreal(Z) & ismatrix(Z), 'stats:classreg:regr:lmeutils:StandardLinearLikeMixedModel:BadZ',num2str(N));
            slme.assertThat(        all(size(Z) == [N,q])          , 'stats:classreg:regr:lmeutils:StandardLinearLikeMixedModel:BadZ',num2str(N));

            % (7) Z must not have any NaN or Inf values.
            slme.assertThat(~slme.hasNaNInf(Z),'stats:classreg:regr:lmeutils:StandardLinearLikeMixedModel:NoNaNInfAllowed','Z');   
            
            % (8) Ensure that Z is sparse.
            if ~issparse(Z)
                Z = sparse(Z);
            end                                 
            
            % (9) Psi is an object of type CovarianceMatrix.
            slme.assertThat(isa(Psi,'classreg.regr.lmeutils.covmats.CovarianceMatrix'),'stats:classreg:regr:lmeutils:StandardLinearLikeMixedModel:BadPsi');
            
            % (10) Psi.Size is q.
            slme.assertThat(Psi.Size == q,'stats:classreg:regr:lmeutils:StandardLinearLikeMixedModel:BadPsi_Size',num2str(q));
            
            % (11) Use subclass specific validation.
            FitMethod = validateFitMethod(slme,FitMethod);            
            
            % (12) N is >= 2. p or q may be 0.
            slme.assertThat(N >= 2,'stats:classreg:regr:lmeutils:StandardLinearLikeMixedModel:BadN');
            
            % (13) N >= p.
            slme.assertThat(N >= p,'stats:classreg:regr:lmeutils:StandardLinearLikeMixedModel:BadN_MinN',num2str(p));             
            
        end % end of validateInputs.        
        
        function slme = invalidateFit(slme)
%invalidateFit - Invalidate the possibly fitted StandardLinearLikeMixedModel.
%   slme = invalidateFit(slme) marks the object as not fitted to data and
%   not ready for stats and returns the modified object.

            % (1) Invalidate the fit.
            slme.isFitToData = false;
            slme.isReadyForStats = false;
            
        end % end of invalidateFit.
        
        function X = validateX(slme,X)
        
            % (1) Ensure X is a numeric real matrix of size slme.N by 
            % slme.p. 
            slme.assertThat(isnumeric(X) & isreal(X) & ismatrix(X),'stats:classreg:regr:lmeutils:StandardLinearLikeMixedModel:InValidX_Size',num2str(slme.N),num2str(slme.p)); 
            slme.assertThat(    all(size(X) == [slme.N,slme.p])   ,'stats:classreg:regr:lmeutils:StandardLinearLikeMixedModel:InValidX_Size',num2str(slme.N),num2str(slme.p)); 
            
            % (2) Ensure no NaN/Inf values in X.
            slme.assertThat(~slme.hasNaNInf(X),'stats:classreg:regr:lmeutils:StandardLinearLikeMixedModel:NoNaNInfAllowed','X');
            
            % (3) Ensure that rank of X is equal to slme.p.
            slme.assertThat(rank(X) == slme.p,'stats:classreg:regr:lmeutils:StandardLinearLikeMixedModel:InValidX_Rank',num2str(slme.p));                        
            
        end % end of validateX.
        
        function y = validatey(slme,y)
            
            % (1) Ensure y is a numeric, real, vector of size slme.N by 1.
            slme.assertThat(isnumeric(y) & isreal(y) & isvector(y),'stats:classreg:regr:lmeutils:StandardLinearLikeMixedModel:InValidY_Size',num2str(slme.N),num2str(1)); 
            slme.assertThat(     all(size(y) == [slme.N,1])       ,'stats:classreg:regr:lmeutils:StandardLinearLikeMixedModel:InValidY_Size',num2str(slme.N),num2str(1)); 

            % (2) Ensure no NaN or Inf values in y.
            slme.assertThat(~slme.hasNaNInf(y),'stats:classreg:regr:lmeutils:StandardLinearLikeMixedModel:NoNaNInfAllowed','y');
            
        end % end of validatey.
        
        function Z = validateZ(slme,Z)
            
            % (1) Ensure Z is a numeric, real matrix of size slme.N by
            % slme.q.
            slme.assertThat(isnumeric(Z) & isreal(Z) & ismatrix(Z),'stats:classreg:regr:lmeutils:StandardLinearLikeMixedModel:InValidZ_Size',num2str(slme.N),num2str(slme.q));
            slme.assertThat(    all(size(Z) == [slme.N,slme.q])   ,'stats:classreg:regr:lmeutils:StandardLinearLikeMixedModel:InValidZ_Size',num2str(slme.N),num2str(slme.q));

            % (2) Ensure no NaN/Inf values in Z.
            slme.assertThat(~slme.hasNaNInf(Z),'stats:classreg:regr:lmeutils:StandardLinearLikeMixedModel:NoNaNInfAllowed','Z');
            
            % (3) Ensure that Z is sparse.
            if ~issparse(Z)
                Z = sparse(Z);
            end                       
            
        end % end of validateZ.
        
        function Psi = validatePsi(slme,Psi)
           
            % (1) Ensure Psi is an object of type CovarianceMatrix.
            slme.assertThat(isa(Psi,'classreg.regr.lmeutils.covmats.CovarianceMatrix'),'stats:classreg:regr:lmeutils:StandardLinearLikeMixedModel:BadPsi');
            
            % (2) Ensure Psi.Size is slme.q.
            slme.assertThat(Psi.Size == slme.q,'stats:classreg:regr:lmeutils:StandardLinearLikeMixedModel:BadPsi_Size',num2str(slme.q));            
            
        end % end of validatePsi.                
        
        function initializationmethod = validateInitializationMethod(slme,initializationmethod)
            
            initializationmethod = internal.stats.getParamVal(initializationmethod,slme.AllowedInitializationMethods,'InitializationMethod');
            
        end % end of validateInitializationMethod.
        
        function theta = validateTheta(slme,theta)

            % (1) Ensure theta is a numeric, real vector.
            assert( isnumeric(theta) & isreal(theta) & isvector(theta) );
            
            % (2) Ensure that theta has the right size.
            assert( all(size(theta) == [slme.Psi.NumParametersExcludingSigma,1]) );
            
            % (3) Ensure no NaN or Inf values in theta.
            assert( ~slme.hasNaNInf(theta) );
            
        end % end of validateTheta.
        
        function beta = validateBeta(slme,beta)

            % (1) Ensure that beta is a numeric, real vector.
            assert( isnumeric(beta) & isreal(beta) & isvector(beta) );
            
            % (2) Ensure that beta has the right size.
            assert( all(size(beta) == [slme.p,1]) );
            
            % (3) Ensure no NaN/Inf values in beta.
            assert( ~slme.hasNaNInf(beta) );
            
        end % end of validateBeta.
        
        function sigma = validateSigma(slme,sigma)

            % (1) Ensure that sigma is a numeric, real scalar.
            assert( isnumeric(sigma) & isreal(sigma) & isscalar(sigma) );
            
            % (2) Ensure that sigma is positive.
            assert( sigma > 0 );
            
            % (3) Ensure no NaN/Inf values in sigma.
            assert( ~slme.hasNaNInf(sigma) );
            
        end % end of validateSigma.
                
        function [optimizer,optimizeroptions,initializationmethod] ...
                    = validateOptimizationOptions(slme,optimizer,optimizeroptions,initializationmethod)
                
            % (1) optimizer must be a string from the list AllowedOptimizers.
            optimizer = internal.stats.getParamVal(optimizer,slme.AllowedOptimizers,'Optimizer');                     
            
            % (2) initializationmethod must be a string from the list AllowedInitializationMethods. 
            initializationmethod = internal.stats.getParamVal(initializationmethod,slme.AllowedInitializationMethods,'InitializationMethod');
                            
            % (3) optimizeroptions must be a struct or a optim.options.SolverOptions object.
            slme.assertThat(isstruct(optimizeroptions) || isa(optimizeroptions,'optim.options.SolverOptions'),'stats:classreg:regr:lmeutils:StandardGeneralizedLinearMixedModel:BadOptimizerOptions');            
            switch lower(optimizer)
                case 'quasinewton'
                    % (1) optimizeroptions can be empty or a struct.
                    slme.assertThat(isempty(optimizeroptions) || isstruct(optimizeroptions),'stats:classreg:regr:lmeutils:StandardGeneralizedLinearMixedModel:BadOptimizerOptions_quasinewton');                                                            
                case {'fminunc'}
                    % (1) optimizeroptions can be empty or of class optim.options.SolverOptions.
                    slme.assertThat(isempty(optimizeroptions) || isa(optimizeroptions,'optim.options.SolverOptions'),'stats:classreg:regr:lmeutils:StandardGeneralizedLinearMixedModel:BadOptimizerOptions_fminunc');  
                case {'fminsearch'}
                    % (1) optimizeroptions can be empty or a struct.
                    slme.assertThat(isempty(optimizeroptions) || isstruct(optimizeroptions),'stats:classreg:regr:lmeutils:StandardGeneralizedLinearMixedModel:BadOptimizerOptions_fminsearch');                    
            end
            
        end % end of validateOptimizationOptions.                
    end
        
    % Hypothesis tests on fixed/random effects.
    methods (Access=public)        
        % Public method for testing c'*beta = e. (T-test).
        function [P,T,DF] = betaTTest(slme,c,e,dfmethod)
%betaTTest - T-test for testing c'*beta = e.            
%   [P,T,DF] = betaTTest(slme,c,e,dfmethod) takes an object of type
%   StandardLinearLikeMixedModel slme, a p by 1 vector c, a scalar e and a
%   string dfmethod which is one of 'none', 'residual' or 'satterthwaite'
%   and performs the T-test for testing H0: c'*beta = e. P is the two sided
%   p-value for the test, T is the test statistic and DF is the denominator
%   degrees of freedom used in the T-test.
    
            % TODO: Account for an all zero c.
            % (1) Ensure that c is sensible.
            if size(c,1) == 1
                c = c';
            end
            assert( all(size(c) == [slme.p,1]) );
            
            % (2) Ensure that e is sensible.            
            assert( all(size(e) == [1,1]) );
            
            % (3) Ensure that dfmethod is one of 'none','residual' or
            % 'satterthwaite'.
            dfmethod = lower(dfmethod);
            assert( ischar(dfmethod) & isvector(dfmethod) ...
                & (size(dfmethod,1) == 1) );
            assert( any(strcmpi(dfmethod,slme.AllowedDFMethods)) );
            
            % (4) Compute the df based on dfmethod.
            switch dfmethod
                case 'none'
                    DF = Inf;
                case 'residual'
                    DF = slme.N - slme.rankX;
                case 'satterthwaite'
                    DF = dfBetaTTest(slme,c);
                otherwise
                    % <entry key="BadDFMethod">Degrees of freedom method must be ''None'', ''Residual'' or ''Satterthwaite''.</entry>
                    error(message('stats:classreg:regr:lmeutils:StandardLinearLikeMixedModel:BadDFMethod'));
            end
            
            % (5) Compute the test statistic T.
            T = (c'*slme.betaHat - e)/sqrt(c'*slme.covbetaHat*c);
            
            % (6) Compute the two sided p-value: P = 2*Pr(t_df > abs(T)).
            % P = 2*(1 - tcdf(abs(T),DF));
            P = 2*(tcdf(abs(T),DF,'upper'));
            
        end % end of betaTTest.
        
        % Public method for testing t'*beta + s'*b = e (T-test).
        function [P,T,DF] = betaBTTest(slme,t,s,e,dfmethod)
%betaBTTest - T-test for testing t'*beta + s'*b = e.            
%   [P,T,DF] = betaBTTest(slme,t,s,e,dfmethod) takes an object of type
%   StandardLinearLikeMixedModel slme, a p by 1 vector t, a q by 1 vector
%   s, a scalar e and a string dfmethod which is one of 'none', 'residual'
%   or 'satterthwaite' and performs the T-test for testing H0: t'*beta +
%   s'*b = e. P is the two sided p-value for the test, T is the test
%   statistic and DF is the denominator degrees of freedom used in the
%   T-test.
            
            % (1) Ensure that t is sensible.
            if size(t,1) == 1
                t = t';
            end
            assert( all(size(t) == [slme.p,1]) );
            
            % (2) Ensure that s is sensible.
            if size(s,1) == 1
                s = s';
            end
            assert( all(size(s) == [slme.q,1]) );
            
            % (3) Ensure that e is sensible.
            assert( all(size(e) == [1,1]) );
            
            % (4) Get the joint contast vector c.
            c = [t;s];                        
                        
            % (5) Ensure that dfmethod is one of 'none','residual' or
            % 'satterthwaite'.
            dfmethod = lower(dfmethod);
            assert( ischar(dfmethod) & isvector(dfmethod) ...
                & (size(dfmethod,1) == 1) );
            assert( any(strcmpi(dfmethod,slme.AllowedDFMethods)) );

            % (6) Compute the df based on dfmethod.
            switch dfmethod
                case 'none'
                    DF = Inf;
                case 'residual'
                    DF = slme.N - slme.rankX;
                case 'satterthwaite'
                    DF = dfBetaBTTest(slme,c);
                otherwise
                    % <entry key="BadDFMethod">Degrees of freedom method must be ''None'', ''Residual'' or ''Satterthwaite''.</entry>
                    error(message('stats:classreg:regr:lmeutils:StandardLinearLikeMixedModel:BadDFMethod'));
            end
            
            % (7) Get the joint vector [betaHat;bHat].
            betaHatbHat = [slme.betaHat;slme.bHat];
                       
            % (8) Compute the test statistic T. We are effectively doing
            % the same thing as: 
            %   V = covBetaHatBHat(slme); 
            %   T = (c'*betaHatbHat - e)/sqrt(c'*V*c);
            % but more efficiently.
            Xnew = t';
            Znew = s';
            [~,ctVc] = getEstimateAndVariance(slme,Xnew,Znew,slme.betaHat,...
                slme.DeltabHat,slme.thetaHat,slme.sigmaHat);
            T = (c'*betaHatbHat - e)/sqrt(ctVc);
                        
            % (9) Compute the two sided p-value: P = 2*Pr(t_df > abs(T)).
            % P = 2*(1 - tcdf(abs(T),DF));   
            P = 2*(tcdf(abs(T),DF,'upper'));

        end % end of betaBTTest.
        
        % Public method for testing L*beta = e. (F-test).
        function [P,T,DF1,DF2] = betaFTest(slme,L,e,dfmethod)
%betaFTest - F-test for testing L*beta = e.            
%   [P,T,DF1,DF2] = betaFTest(slme,L,e,dfmethod) takes an object of type
%   StandardLinearLikeMixedModel slme, a r by p matrix L of rank r, a r by
%   1 vector e and a string dfmethod which is one of 'none', 'residual' or
%   'satterthwaite' and performs the F-test for testing H0: L*beta = e. P
%   is the p-value for the test, T is the test statistic, DF1 is the
%   numerator degrees of freedom and DF2 is the denominator degrees of
%   freedom used in the F-test.

            % (1) Ensure that L is sensible.
            assert( size(L,2) == slme.p );

            % (2) Ensure that e is sensible.
            if size(e,1) == 1
                e = e';
            end
            assert( all(size(e) == [size(L,1),1]) );
            
            % (3) Handle the empty L or e case.
            if isempty(L)
                  P = NaN;
                  T = NaN;
                DF1 = NaN;
                DF2 = NaN;
                return;
            end
            
            % (4) Get effective L and e.
            %
            % L may not be of full row rank and the hypothesis L*beta = e 
            % may not be consistent. Get the effective L and e for the test
            % and ensure that original L and e are consistent.            
            %       L = r by slme.p matrix
            %       e = r by 1 vector.
            % Desired test: 
            %       L*beta = e.            
            [ok,~,L,e] = slme.fullrankH(L,e);                        
            if ~ok
                % <entry key="BadHypothesisMatrix">The hypothesis matrix and test value are not consistent.</entry>
                error(message('stats:classreg:regr:lmeutils:StandardLinearLikeMixedModel:BadHypothesisMatrix'));
            end
            r = size(L,1);
            
            % (5) Ensure that dfmethod is one of 'none','residual' or
            % 'satterthwaite'.
            dfmethod = lower(dfmethod);
            assert( ischar(dfmethod) & isvector(dfmethod) ...
                & (size(dfmethod,1) == 1) );
            assert( any(strcmpi(dfmethod,slme.AllowedDFMethods)) );
            
            % (6) Compute the denominator df based on dfmethod.
            switch dfmethod
                case 'none'
                    DF2 = Inf;
                case 'residual'
                    DF2 = slme.N - slme.rankX;
                case 'satterthwaite'
                    DF2 = dfBetaFTest(slme,L);
                otherwise
                    % <entry key="BadDFMethod">Degrees of freedom method must be ''None'', ''Residual'' or ''Satterthwaite''.</entry>
                    error(message('stats:classreg:regr:lmeutils:StandardLinearLikeMixedModel:BadDFMethod'));
            end
            
            % (7) Compute the test statistic T.
            delta = (L*slme.betaHat - e);
            T = delta'*( (L*slme.covbetaHat*L') \ delta )/r;
            
            % (8) Set numerator degrees of freedom as r.
            DF1 = r;
            
            % (9) Compute p-value for the F-test.
            % P = 1 - fcdf(T,DF1,DF2);
            P = fcdf(T,DF1,DF2,'upper');
            
        end % end of betaFTest.

        % Public method for testing H*beta + K*b = e. (F-test).
        function [P,T,DF1,DF2] = betaBFTest(slme,H,K,e,dfmethod)
%betaBFTest - F-test for testing H*beta + K*b = e.            
%   [P,T,DF1,DF2] = betaBFTest(slme,H,K,e,dfmethod) takes an object of type
%   StandardLinearLikeMixedModel slme, a r by p matrix H, a r by q matrix
%   K, a r by 1 vector e and a string dfmethod which is one of 'none',
%   'residual' or 'satterthwaite' and performs the F-test for testing H0:
%   H*beta + K*b = e. P is the p-value for the test, T is the test
%   statistic, DF1 is the numerator degrees of freedom and DF2 is the
%   denominator degrees of freedom used in the F-test. The rank of matrix
%   [H,K] must be r <= (p + q).

            % (1) Ensure that H is sensible.
            assert( size(H,2) == slme.p );
                        
            % (2) Ensure that K is sensible.
            assert( size(K,1) == size(H,1) );
            assert( size(K,2) == slme.q );
            
            % (3) Form contrast matrix L for [beta;b].
            L = [H,K];                        
            
            % (4) Ensure that e is sensible.
            if size(e,1) == 1
                e = e';
            end
            assert( all(size(e) == [size(H,1),1]) );
            
            % (5) Get effective L and e.
            %
            % L may not be of full row rank and the hypothesis L*[beta;b] =
            % e may not be consistent. Get the effective L and e for the
            % test and ensure that original L and e are consistent.
            %       L = r by [slme.p + slme.q] matrix
            %       e = r by 1 vector.
            % Desired test: 
            %       L*[beta;b] = e.            
            [ok,~,L,e] = slme.fullrankH(L,e);                        
            if ~ok
                % <entry key="BadHypothesisMatrix">The hypothesis matrix and test value are not consistent.</entry>               
                error(message('stats:classreg:regr:lmeutils:StandardLinearLikeMixedModel:BadHypothesisMatrix'));
            end
            r = size(L,1);            
            
            % (6) Ensure that dfmethod is one of 'none','residual' or
            % 'satterthwaite'.
            dfmethod = lower(dfmethod);
            assert( ischar(dfmethod) & isvector(dfmethod) ...
                & (size(dfmethod,1) == 1) );
            assert( any(strcmpi(dfmethod,slme.AllowedDFMethods)) );
            
            % (7) Compute the denominator df based on dfmethod.
            switch dfmethod
                case 'none'
                    DF2 = Inf;
                case 'residual'
                    DF2 = slme.N - slme.rankX;
                case 'satterthwaite'
                    DF2 = dfBetaBFTest(slme,L);
                otherwise
                    % <entry key="BadDFMethod">Degrees of freedom method must be ''None'', ''Residual'' or ''Satterthwaite''.</entry>
                    error(message('stats:classreg:regr:lmeutils:StandardLinearLikeMixedModel:BadDFMethod'));
            end
            
            % (8) Get the joint vector [betaHat;bHat].
            betaHatbHat = [slme.betaHat;slme.bHat];                        
                                    
            % (9) Compute the test statistic T. We are effectively doing
            % the same thing as:
            %   V = covBetaHatBHat(slme);
            %   T = delta'*( (L*V*L') \ delta )/r;
            % but more efficiently.
            delta = (L*betaHatbHat - e);            
            p = slme.p;
            q = slme.q;
            Xnew = L(:,1:p);
            Znew = L(:,p+1:p+q);
            [~,LVLt] = getEstimateAndVariance(slme,Xnew,Znew,slme.betaHat,...
                slme.DeltabHat,slme.thetaHat,slme.sigmaHat,'covariance');
            T = delta'*( LVLt \ delta )/r;
            
            % (10) Set numerator degrees of freedom as r.
            DF1 = r;
            
            % (11) Compute p-value for the F-test.
            % P = 1 - fcdf(T,DF1,DF2);
            P = fcdf(T,DF1,DF2,'upper');
            
        end % end of betaBFTest.
    end            
    
    % Confidence intervals on fixed/random effects.
    methods (Access=public)
        % Public method for computing CI on c'*beta.
        function [CI,DF] = betaCI(slme,c,alpha,dfmethod)
%betaCI - Confidence interval on c'*beta.            
%   [CI,DF] = betaCI(slme,c,alpha,dfmethod) takes an object of type
%   StandardLinearLikeMixedModel slme, a p by 1 vector c, a scalar alpha
%   between 0 and 1 and a string dfmethod which is one of 'none',
%   'residual' or 'satterthwaite' and computes a (1-alpha) confidence
%   interval for c'*beta. CI has two columns and 1 row. CI(1) and CI(2)
%   contain respectively the lower and upper limit of the confidence
%   interval for c'*beta. DF is the degrees of freedom used to compute the
%   confidence interval.
   
            % (1) Ensure that c is sensible.
            if size(c,1) == 1
                c = c';
            end
            assert( all(size(c) == [slme.p,1]) );
            
            % (2) Ensure that alpha is sensible.
            assert( isnumeric(alpha) & isreal(alpha) & isscalar(alpha) );
            assert( alpha >= 0 & alpha <= 1 );
                        
            % (3) Ensure that dfmethod is one of 'none','residual' or
            % 'satterthwaite'.
            dfmethod = lower(dfmethod);
            assert( ischar(dfmethod) & isvector(dfmethod) ...
                & (size(dfmethod,1) == 1) );
            assert( any(strcmpi(dfmethod,slme.AllowedDFMethods)) );
            
            % (4) Compute the df based on dfmethod.
            switch dfmethod
                case 'none'
                    DF = Inf;
                case 'residual'
                    DF = slme.N - slme.rankX;
                case 'satterthwaite'
                    DF = dfBetaTTest(slme,c);
                otherwise
                    % <entry key="BadDFMethod">Degrees of freedom method must be ''None'', ''Residual'' or ''Satterthwaite''.</entry>
                    error(message('stats:classreg:regr:lmeutils:StandardLinearLikeMixedModel:BadDFMethod'));                    
            end
            
            % (5) Get the covariance of betaHat.
            V = slme.covbetaHat;
            
            % (6) Compute half width of the CI.
            delta = tinv(1-alpha/2,DF)*sqrt(c'*V*c);
            
            % (7) Form the output CI.
            ctbetaHat = c'*slme.betaHat;
            CI = [ctbetaHat - delta, ctbetaHat + delta];
            
        end % end of betaCI.
        
        % Public method for computing CI on t'*beta + s'*b.
        function [CI,DF] = betaBCI(slme,t,s,alpha,dfmethod)
%betaBCI - Confidence interval on t'*beta + s'*b.            
%   [CI,DF] = betaBCI(slme,t,s,alpha,dfmethod) takes an object of type
%   StandardLinearLikeMixedModel slme, a p by 1 vector t, a q by 1 vector
%   s, a scalar alpha between 0 and 1 and a string dfmethod which is one of
%   'none', 'residual' or 'satterthwaite' and computes a (1-alpha)
%   confidence interval for t'*beta + s'*b. CI has two columns and 1 row.
%   CI(1) and CI(2) contain respectively the lower and upper limit of the
%   confidence interval for t'*beta + s'*b. DF is the degrees of freedom
%   used to compute the confidence interval.
            
            % (1) Ensure that t is sensible.
            if size(t,1) == 1
                t = t';
            end
            assert( all(size(t) == [slme.p,1]) );
            
            % (2) Ensure that s is sensible.
            if size(s,1) == 1
                s = s';
            end
            assert( all(size(s) == [slme.q,1]) );
            
            % (3) Get the joint contast vector c.
            c = [t;s];
            
            % (3) Ensure that alpha is sensible.
            assert( isnumeric(alpha) & isreal(alpha) & isscalar(alpha) );
            assert( alpha >= 0 & alpha <= 1 );
                        
            % (4) Ensure that dfmethod is one of 'none','residual' or
            % 'satterthwaite'.
            dfmethod = lower(dfmethod);
            assert( ischar(dfmethod) & isvector(dfmethod) ...
                & (size(dfmethod,1) == 1) );
            assert( any(strcmpi(dfmethod,slme.AllowedDFMethods)) );

            % (5) Compute the df based on dfmethod.
            switch dfmethod
                case 'none'
                    DF = Inf;
                case 'residual'
                    DF = slme.N - slme.rankX;
                case 'satterthwaite'
                    DF = dfBetaBTTest(slme,c);
                otherwise
                    % <entry key="BadDFMethod">Degrees of freedom method must be ''None'', ''Residual'' or ''Satterthwaite''.</entry>
                    error(message('stats:classreg:regr:lmeutils:StandardLinearLikeMixedModel:BadDFMethod'));                    
            end                                   

            % (6) Compute half width of the CI. We are effectively doing
            % the same thing as:
            %   V = covBetaHatBHat(slme);
            %   delta = tinv(1-alpha/2,DF)*sqrt(c'*V*c);
            % but more efficiently.
            Xnew = t';
            Znew = s';
            [~,ctVc] = getEstimateAndVariance(slme,Xnew,Znew,slme.betaHat,...
                slme.DeltabHat,slme.thetaHat,slme.sigmaHat);            
            delta = tinv(1-alpha/2,DF)*sqrt(ctVc);
            
            % (7) Form the output CI.
            ctbetaHatbHat = c'*[slme.betaHat;slme.bHat];
            CI = [ctbetaHatbHat - delta, ctbetaHatbHat + delta];
            
        end % end of betaBCI.                
    end
    
    % Table of fixed/random effects and predictions.
    methods (Access=public)    
        function fetable = fixedEffects(slme,alpha,dfmethod)
%fixedEffects - Make a table of fixed effects.
%   fetable = fixedEffects(slme,alpha,dfmethod) takes an object slme of
%   type StandardLinearLikeMixedModel and makes a fixed effects table.
%   fetable is a table array with the following columns:
%
%           Estimate - Estimate of the fixed effect
%           SE - Standard deviation of the estimate of fixed effect.
%           tStat - t-stat for testing the true fixed effect against 0.
%           DF - Degrees of freedom used in the t-test. 
%           pValue - p-value from the t-test.
%           Lower - Lower limit of a (1-alpha) CI.
%           Upper - Upper limit of a (1-alpha) CI.
%
%   The input alpha is a number between 0 and 1 and dfmethod is one of
%   'satterthwaite','none' or 'residual'.

            % (1) Ensure that alpha is sensible.
            assert( isnumeric(alpha) & isreal(alpha) & isscalar(alpha) );
            assert( alpha >= 0 & alpha <= 1 );
                        
            % (2) Ensure that dfmethod is one of 'none','residual' or
            % 'satterthwaite'.
            dfmethod = lower(dfmethod);
            assert( ischar(dfmethod) & isvector(dfmethod) ...
                & (size(dfmethod,1) == 1) );
            assert( any(strcmpi(dfmethod,slme.AllowedDFMethods)) );

            % (3) How many fixed effects?
            p = slme.p;
            
            % (4) Get predicted beta and its se.
            pred = slme.betaHat;
            if (p == 0)
                % No fixed effects case.
                se = zeros(p,1);
            else
                se = sqrt(diag(slme.covbetaHat));
            end
            
            % (5) Compute denominator degrees of freedom.
            switch dfmethod
                case 'none'
                    DF = Inf*ones(p,1);
                case 'residual'                    
                    DF = (slme.N - p)*ones(p,1);
                case 'satterthwaite'
                    DF = zeros(p,1);
                    Ip = speye(p);
                    for i = 1:p
                        % Make a p-by-1 contrast vector.
                        c = full(Ip(:,i));
                        DF(i) = dfBetaTTest(slme,c);
                    end
                otherwise                    
                    error(message('stats:classreg:regr:lmeutils:StandardLinearLikeMixedModel:BadDFMethod'));
            end
            
            % (6) Compute the test statistic T.
            T = pred ./ se;
            
            % (7) Compute the two sided p-value: P = 2*Pr(t_df > abs(T)).
            % P = 2*(1 - tcdf(abs(T),DF));   
            P = 2*(tcdf(abs(T),DF,'upper'));
                        
            % (8) Compute half width of the CI.
            halfwidth = tinv(1-alpha/2,DF) .* se;
            
            % (9) Form the output CI.        
            LB = pred - halfwidth;
            UB = pred + halfwidth;                       
            
            % (10) Make sure that none of the DF values are NaN.
            if any(isnan(DF))
                error(message('stats:classreg:regr:lmeutils:StandardLinearLikeMixedModel:ErrorDFCalculation'));                
            end                                    
            
            % (11) Make the output table array.
            fetable = table(pred,se,T,DF,P,LB,UB,'VariableNames',{'Estimate','SE',...
                'tStat','DF','pValue','Lower','Upper'});
            
        end % end of fixedEffects.            
        
        function retable = randomEffects(slme,alpha,dfmethod)
%randomEffects - Make a table of random effects.
%   retable = randomEffects(slme,alpha,dfmethod) takes an object slme of
%   type StandardLinearLikeMixedModel and makes a random effects table.
%   retable is a table array with the following columns:
%
%           Estimate - Estimate of the random effect
%           SEPred - Standard deviation of elements of (bHat - b).
%           tStat - t-stat for testing the random effect against 0.
%           DF - Degrees of freedom used in the t-test. 
%           pValue - p-value from the t-test.
%           Lower - Lower limit of a (1-alpha) CI.
%           Upper - Upper limit of a (1-alpha) CI.
%
%   The input alpha is a number between 0 and 1 and dfmethod is one of
%   'satterthwaite','none' or 'residual'.            
            
             % (1) Ensure that alpha is sensible.
            assert( isnumeric(alpha) & isreal(alpha) & isscalar(alpha) );
            assert( alpha >= 0 & alpha <= 1 );
                        
            % (2) Ensure that dfmethod is one of 'none','residual' or
            % 'satterthwaite'.
            dfmethod = lower(dfmethod);
            assert( ischar(dfmethod) & isvector(dfmethod) ...
                & (size(dfmethod,1) == 1) );
            assert( any(strcmpi(dfmethod,slme.AllowedDFMethods)) );

            % (3) Get p and q from slme.
            p = slme.p;
            q = slme.q;                       
            
            % (4) Get the SE Pred column. Create Xnew = zeros(q,p) and 
            % Znew = eye(q) to extract elements of bHat as pred and its
            % estimated variance in varpred.
            Xnew = sparse(q,p);
            Znew = speye(q);
            [pred,varpred] = getEstimateAndVariance(slme,Xnew,Znew,slme.betaHat,slme.DeltabHat,slme.thetaHat,slme.sigmaHat);
            se = sqrt(varpred);
            
            % (5) Compute denominator degrees of freedom.
            switch dfmethod
                case 'none'
                    DF = Inf*ones(q,1);
                case 'residual'
                    DF = (slme.N - slme.p)*ones(q,1);
                case 'satterthwaite'
                    DF = zeros(q,1);
                    for i = 1:q
                        % Make a (p+q)-by-1 contrast vector.
                        c = full([Xnew(i,:),Znew(i,:)])';
                        DF(i) = dfBetaBTTest(slme,c);
                    end
                otherwise                    
                    error(message('stats:classreg:regr:lmeutils:StandardLinearLikeMixedModel:BadDFMethod'));
            end
            
            % (6) Compute the test statistic T.
            T = pred ./ se;
            
            % (7) Compute the two sided p-value: P = 2*Pr(t_df > abs(T)).
            % P = 2*(1 - tcdf(abs(T),DF));   
            P = 2*(tcdf(abs(T),DF,'upper'));
                        
            % (8) Compute half width of the CI.
            halfwidth = tinv(1-alpha/2,DF) .* se;
            
            % (9) Form the output CI.        
            LB = pred - halfwidth;
            UB = pred + halfwidth;                       
            
            % (10) Make sure that none of the DF values are NaN.
            if any(isnan(DF))
                error(message('stats:classreg:regr:lmeutils:StandardLinearLikeMixedModel:ErrorDFCalculation'));                
            end            
            
            % (11) Make the output table array.
            retable = table(pred,se,T,DF,P,LB,UB,'VariableNames',{'Estimate','SEPred',...
                'tStat','DF','pValue','Lower','Upper'});
            
        end % end of randomEffects.
                        
        function predtable = predictTable(slme,Xnew,Znew,alpha,dfmethod)
%predictTable - Predict output from a fitted StandardLinearLikeMixedModel.
%   predtable = predictTable(slme,Xnew,Znew,alpha,dfmethod) takes a fitted
%   StandardLinearLikeMixedModel slme, a fixed effects design matrix Xnew,
%   a random effects design matrix Znew, a confidence level alpha, a
%   degrees of freedom method dfmethod and produces (1-alpha) pointwise
%   confidence intervals for the predictions from slme at Xnew and Znew.
%   The predictions include contributions from both fixed effects and BLUPs
%   of random effects. Xnew must be M by p and Znew must be M by q where p
%   = slme.p and q = slme.q. The output predtable is a table array with the
%   following columns:
%
%       Pred - The prediction
%       SEPred - SE of Prediction error
%       DF - The degrees of freedom used to compute the confidence interval
%       Alpha - The alpha value used for the confidence interval
%       Lower - Lower limit of the confidence interval
%       Upper - Upper limit of the confidence interval

            % (1) Ensure that Xnew is sensible.
            assert( isnumeric(Xnew) & isreal(Xnew) & ismatrix(Xnew) );
            assert( size(Xnew,2) == slme.p );
            
            % (2) Ensure that Znew is sensible.
            assert( isnumeric(Znew) & isreal(Znew) & ismatrix(Znew) );
            assert( size(Znew,2) == slme.q );
            
            % (3) Ensure that Xnew and Znew have the same number of rows.
            assert( size(Xnew,1) == size(Znew,1) );

            % (4) Ensure that alpha is sensible.
            assert( isnumeric(alpha) & isreal(alpha) & isscalar(alpha) );
            assert( alpha >= 0 & alpha <= 1 );
                        
            % (5) Ensure that dfmethod is one of 'none','residual' or
            % 'satterthwaite'.
            dfmethod = lower(dfmethod);
            assert( ischar(dfmethod) & isvector(dfmethod) ...
                & (size(dfmethod,1) == 1) );
            assert( any(strcmpi(dfmethod,slme.AllowedDFMethods)) );
      
            % (6) Get the number of rows M in Xnew.            
            M = size(Xnew,1);
                       
            % (7) Get the covariance of [betaHat-beta;bHat-b].
            V = covBetaHatBHat(slme);                        
            
            % (8) Get betaHat and bHat.
            betaHat = slme.betaHat;
            bHat = slme.bHat;
            
            % (9) Get ypred.
            ypred = Xnew*betaHat + Znew*bHat;
            
            % (10) Get Alpha.
            Alpha = alpha*ones(M,1);
            
            % (11) Initialize Pred, SEPred, DF, Alpha, LB and UB.
            SEPred = zeros(M,1);
            DF = zeros(M,1);
            LB = zeros(M,1);
            UB = zeros(M,1);                       
            for i = 1:M
                % (12) Compute (1-alpha) CI for t'*beta + s'*b where t is 
                % given by the rows of Xnew and s is given by the rows of 
                % Znew. Store the results in LB, UB and DF.
                t = Xnew(i,:)';
                s = Znew(i,:)';
                [CI0,DF0] = betaBCI(slme,t,s,alpha,dfmethod);
                LB(i) = CI0(1);
                UB(i) = CI0(2);
                DF(i) = DF0;
                
                % (13) Compute sqrt(c'*V*c) where c = [t;s].
                c = [t;s];
                SEPred(i) = sqrt(c'*V*c);                
            end
            
            % (14) Make sure that none of the DF values are NaN.
            if any(isnan(DF))
                % <entry key="ErrorDFCalculation">Degrees of freedom calculation produced a NaN result. Try using a different degrees of freedom method.</entry>
                error(message('stats:classreg:regr:lmeutils:StandardLinearLikeMixedModel:ErrorDFCalculation'));                
            end            
            
            % (15) Make the output table array.
            predtable = table(ypred,SEPred,DF,Alpha,LB,UB,'VariableNames',{'Pred','SEPred','DF','Alpha','Lower','Upper'});   
        end % end of predictTable.        
    
        function [ypred,CI,DF] = predict(slme,Xnew,Znew,alpha,dfmethod,wantConditional,wantPointwise,wantCurve,hasIntercept)
%predict - Predictions for linear predictor of a fitted StandardLinearLikeMixedModel.
%   [ypred,CI,DF] = predict(slme,Xnew,Znew,alpha,dfmethod,wantConditional,wantPointwise,wantCurve,hasIntercept)
%   takes a fitted StandardLinearLikeMixedModel slme, a fixed effects
%   design matrix Xnew, a random effects design matrix Znew, a confidence
%   level alpha, a degrees of freedom method dfmethod, flags for
%   conditional, pointwise and curve predictions and outputs predictions
%   and (1-alpha) confidence intervals (CIs) for the true linear predictor
%   of slme at Xnew and Znew. If wantConditional is true then the
%   predictions include contributions from both fixed effects and estimates
%   of random effects otherwise the predictions include contributions from
%   only the fixed effects. If wantPointwise is true then pointwise CIs are
%   computed otherwise simultaneous CIs are computed. If wantCurve is true
%   then CIs are for the curve excluding variability due to observation
%   noise otherwise variability due to observation noise is included. Xnew
%   must be M by p and Znew must be M by q where p = slme.p and q = slme.q.
%   The output ypred is a M-by-1 vector containing the predictions
%   corresponding to the rows of Xnew and Znew and output CI is a M-by-2
%   matrix. Each row of CI is a (1-alpha) confidence interval for the true
%   prediction such that CI(:,1) contains the lower confidence limit and
%   CI(:,2) contains the upper confidence limit. DF is a M-by-1 vector
%   containing the DF values used in computing the CIs if wantPointwise is
%   true. If wantPointwise is false then DF is a scalar containing the DF
%   value used for computing the Scheffe simultaneous CIs. hasIntercept
%   should be set to true if the fixed effects part of the model has an
%   intercept and false otherwise.
            
            % (1) Ensure that Xnew is sensible.
            assert( isnumeric(Xnew) & isreal(Xnew) & ismatrix(Xnew) );
            assert( size(Xnew,2) == slme.p );
            
            % (2) Ensure that Znew is sensible.
            assert( isnumeric(Znew) & isreal(Znew) & ismatrix(Znew) );
            assert( size(Znew,2) == slme.q );
            
            % (3) Ensure that Xnew and Znew have the same number of rows.
            assert( size(Xnew,1) == size(Znew,1) );

            % (4) Ensure that alpha is sensible.
            assert( isnumeric(alpha) & isreal(alpha) & isscalar(alpha) );
            assert( alpha >= 0 & alpha <= 1 );
                        
            % (5) Ensure that dfmethod is one of 'none','residual' or
            % 'satterthwaite'.
            dfmethod = lower(dfmethod);
            assert( ischar(dfmethod) & isvector(dfmethod) ...
                & (size(dfmethod,1) == 1) );
            assert( any(strcmpi(dfmethod,slme.AllowedDFMethods)) );

            % (6) Ensure that wantConditional is a scalar logical.
            assert( isscalar(wantConditional) & islogical(wantConditional) );
            
            % (7) Ensure that wantPointwise is a scalar logical.
            assert( isscalar(wantPointwise) & islogical(wantPointwise) );
            
            % (8) Ensure that wantCurve is a scalar logical.
            assert( isscalar(wantCurve) & islogical(wantCurve) );
            
            % (9) Ensure that hasIntercept is a scalar logical.
            assert( isscalar(hasIntercept) & islogical(hasIntercept) );
            
            % (10) Get conditional or marginal prediction ypred.
            if ( wantConditional == true )
                ypred = Xnew*slme.betaHat + Znew*slme.bHat;                                
            else
                ypred = Xnew*slme.betaHat;                
            end
            
            if nargout > 1                
                % (11) Get a vector of estimated variance of predictions. 
                % To do this form matrices C and V such that:
                %
                % (a) For marginal predictions:
                %       C = [Xnew] and V = slme.covbetaHat 
                % (b) For conditional predictions:
                %       C = [Xnew,Znew] and V = covBetaHatBHat(slme);
                % 
                % The ith element of the vector of estimated prediction 
                % variances is then given by: 
                % 
                % varpred(i) = C(i,:)*V*C(i,:)'
                %
                % Add variability due to observation noise to varpred if
                % wantCurve is false.                
                if ( wantConditional == true ) 
                    C = [Xnew,Znew];
                    % This is the same as doing:
                    %   V = covBetaHatBHat(slme);
                    %   varpred = sum(C .* (C*V'),2);
                    % but more efficiently.
                    [~,varpred] = getEstimateAndVariance(slme,Xnew,Znew,slme.betaHat,slme.DeltabHat,slme.thetaHat,slme.sigmaHat);
                else
                    C = Xnew;
                    V = slme.covbetaHat;
                    varpred = sum(C .* (C*V'),2);
                end                
                
                % (12) If wantCurve is false, add residual noise variance 
                % slme.sigmaHat^2 to each element of varpred.
                if ( wantCurve == false )
                    varpred = varpred + (slme.sigmaHat)^2;
                end                
                
                % (13) Get the number of rows in Xnew, the effective number
                % of columns in slme.X and the number of columns in slme.Z.
                M = size(Xnew,1);
                rankX = slme.rankX;
                q = slme.q;
                
                % (14) Compute elementwise multiplier crit of sqrt(varpred) 
                % depending on whether wantPointwise is true or false.
                if ( wantPointwise == true )
                    % (1) Get a vector of DF values.
                    switch dfmethod
                        case 'none'
                            DF = Inf*ones(M,1);
                        case 'residual'
                            DF = (slme.N - rankX)*ones(M,1);
                        case 'satterthwaite'
                            DF = zeros(M,1);
                            for i = 1:M
                                if ( wantConditional == true )
                                    DF(i) = dfBetaBTTest(slme,C(i,:));
                                else
                                    DF(i) = dfBetaTTest(slme,C(i,:));
                                end
                            end
                        otherwise
                            % <entry key="BadDFMethod">Degrees of freedom method must be ''None'', ''Residual'' or ''Satterthwaite''.</entry>                            
                            error(message('stats:classreg:regr:lmeutils:StandardLinearLikeMixedModel:BadDFMethod'));
                    end
                    % (2) Compute multiplier vector crit.
                    crit = tinv(1-alpha/2,DF);
                else
                    % (1) DF will be a scalar.                    
                    switch dfmethod
                        case 'none'
                            DF = Inf;
                        case 'residual'
                            DF = (slme.N - rankX);
                        case 'satterthwaite'
                            if ( wantConditional == true )
                                DF = dfBetaBFTest(slme,eye(rankX + q));
                            else
                                DF = dfBetaFTest(slme,eye(rankX));
                            end
                        otherwise
                            % <entry key="BadDFMethod">Degrees of freedom method must be ''None'', ''Residual'' or ''Satterthwaite''.</entry>                            
                            error(message('stats:classreg:regr:lmeutils:StandardLinearLikeMixedModel:BadDFMethod'));
                    end
                    
                    % (2) Set numerator DF for Scheffe multiplier. Add 1 to
                    % numerator df if wantCurve is false (observation).
                    if ( wantConditional == true )
                        if ( wantCurve == true )
                            numDF = rankX + q;
                        else
                            if hasIntercept
                                numDF = rankX + q;
                            else
                                numDF = rankX + q + 1;
                            end
                        end
                    else
                        if ( wantCurve == true )
                            numDF = rankX;
                        else
                            if hasIntercept
                                numDF = rankX;
                            else
                                numDF = rankX + 1;
                            end
                        end
                    end
                    
                    % (3) Compute crit for simultaneous inference.
                    crit = sqrt( numDF*finv(1-alpha,numDF,DF) );                    
                end
                
                % (15) Compute the CIs.
                delta = crit.*sqrt(varpred);
                CI = [ypred - delta, ypred + delta];            
            end
            
        end % end of predict.    
    end
            
    % Store/unstore the covariance of [betaHat - beta;bHat - b] in object.
    methods (Access=public)        
        function slme = storeCovBetaHatBHat(slme)
%storeCovBetaHatBHat - Store the covariance of [betaHat - beta;bHat - b] in object.
%   slme = storeCovBetaHatBHat(slme) takes an object slme of type
%   StandardLinearLikeMixedModel and stores the covariance of [betaHat -
%   beta;bHat - b] in the object and returns the modified object.

            slme.covbetaHatbHat = covBetaHatBHat(slme);

        end % end of storeCovBetaHatBHat.
        
        function slme = unstoreCovBetaHatBHat(slme)
%unstoreCovBetaHatBHat - Remove the stored covariance of [betaHat - beta;bHat - b] from object.
%   slme = unstoreCovBetaHatBHat(slme) takes an object slme of type
%   StandardLinearLikeMixedModel and removes the stored covariance of
%   [betaHat - beta;bHat - b] from the object and returns the modified
%   object.

            slme.covbetaHatbHat = [];

        end % end of storeCovBetaHatBHat.        
    end
        
    % Get various residual types (convenience function).
    methods (Access=public)        
        function r = residuals(slme,wantConditional,residualType)
%residuals - Convenience function to get residuals of various types.
%   r = residuals(slme,wantConditional,residualType) takes an object slme
%   of type StandardLinearLikeMixedModel and returns the conditional
%   residuals if wantConditional is true and marginal residuals if
%   wantConditional is false. residualType is a string indicating the type
%   of residuals required. residualType can be either 'Raw', 'Pearson',
%   'Standardized'. Note that 'Standardized' will give the internally
%   studentized residuals.

            % (1) Ensure that wantConditional is sensible.
            assert( isscalar(wantConditional) & islogical(wantConditional) );
            
            % (2) Ensure that residualType is one of 'Raw', 'Pearson' or
            % 'Standardized'.           
            residualType = lower(residualType);
            assert( ischar(residualType) & isvector(residualType) ...
                & (size(residualType,1) == 1) );
            assert( any(strcmpi(residualType,slme.AllowedResidualTypes)) );
            
            % (3) Get the desired type of residual.
            switch residualType                
                case 'raw'
                    r = getRawResiduals(slme,wantConditional);
                case 'pearson'
                    r = getPearsonResiduals(slme,wantConditional);
                case 'standardized'
                    % 'standardized' really means internally studentized.
                    % This is what getStudentizedResiduals returns.
                    r = getStudentizedResiduals(slme,wantConditional);
                otherwise
                    error(message('stats:classreg:regr:lmeutils:StandardLinearLikeMixedModel:BadResidualType'));                    
            end
            
        end % end of residuals.        
    end
               
    % Stats on covariance parameters.
    methods (Access=public)                
        % Table of covariance parameters.
        function covtable = covarianceParameters(slme,alpha,wantCIs)
%covarianceParameters - Compute a table containing covariance parameters.
%   covtable = covarianceParameters(slme,alpha,wantCIs) takes an object
%   slme of type StandardLinearLikeMixedModel, alpha which is a number
%   between 0 and 1 and a flag wantCIs that is either true or false and 
%   returns a table array with the following columns:
%
%           Estimate - Estimate of the covariance parameter 
%           Lower - Lower bound of (1-alpha) CI for the covariance parameter
%           Upper - Upper bound of (1-alpha) CI for the covariance parameter
%            
%   The last row of the table corresponds to the residual standard
%   deviation or the square root of the dispersion parameter. If wantCIs is
%   false, the CIs are not returned. You can omit wantCIs in which case
%   wantCIs will be set to true.

            % (0) If called with 2 args, set wantCIs to true.
            if nargin < 3
                wantCIs = true;
            end

            % (1) Ensure that alpha is sensible.
            assert( isnumeric(alpha) & isreal(alpha) & isscalar(alpha) );
            assert( alpha >= 0 & alpha <= 1 );
            
            % (2) Get the estimate of canonical parameters heta.
            hetaHat = slme.Psi.getCanonicalParameters;
            
            % (3) Get the estimate of sigma.
            sigmaHat = slme.sigmaHat;
            
            % (4) Form x = [hetaHat;sigmaHat].
            x = [hetaHat;sigmaHat];

            % (5) Create output covtable with just the 'Estimate'.
            covtable = table(x,'VariableNames',{'Estimate'});
            
            % (6) Get (1-alpha) CI for [heta;sigma]. This may not work
            % sometimes because the covariance matrix of [etaHat;log(sigmaHat)]
            % is not positive definite. When this happens we error out. 
            if wantCIs == true
                try                
                    CI = hetaSigmaCI(slme,alpha);
                    covtable.Lower = CI(:,1);
                    covtable.Upper = CI(:,2);
                catch ME
                    msgID = 'stats:classreg:regr:lmeutils:StandardLinearLikeMixedModel:Message_TooManyCovarianceParameters';
                    msgStr = getString(message(msgID));
                    baseME = MException(msgID,msgStr);
                    ME = addCause(baseME,ME);
                    warning('stats:classreg:regr:lmeutils:StandardLinearLikeMixedModel:Message_TooManyCovarianceParameters',ME.message);
                end                
            end
                        
        end % end of covarianceParameters.        
    
        % CIs on Natural parameter vector.
        function CI = etaLogSigmaCI(slme,alpha)
%etaLogSigmaCI - Confidence intervals on the Natural parameter vector.
%   CI = etaLogSigmaCI(slme,alpha) takes a StandardLinearLikeMixedModel
%   object slme and a confidence level between 0 and 1 and computes the
%   (1-alpha)*100 CIs for [eta;log(sigma)]. CI(:,1) contains the lower
%   limits of the CIs and CI(:,2) contains the upper limits of the CIs. The
%   last row of CI contains the lower and upper limit of the CI for
%   log(sigma).

            % (1) Ensure that alpha is sensible.
            assert( isnumeric(alpha) & isreal(alpha) & isscalar(alpha) );
            assert( alpha >= 0 & alpha <= 1 );
            
            % (2) Form x = [etaHat;log(sigmaHat)].
            etaHat = getNaturalParameters(slme.Psi);
            sigmaHat = slme.sigmaHat;
            x = [etaHat;log(sigmaHat)];
            
            % (3) Get covariance of [etaHat;log(sigmaHat)]. Create the 
            % covariance on the fly if required.           
            C = slme.covetaHatlogsigmaHat;
            if isempty(C)
                C = covEtaHatLogSigmaHat(slme); 
            end
            
            % (4) (1-alpha) confidence interval on elements of 
            % [eta;log(sigma)]. CI will have length(etaHat) + 1 rows and 2 
            % columns.
            % delta = norminv(1-alpha/2)*sqrt(diag(C));
            delta = -norminv(alpha/2)*sqrt(diag(C));
            CI = [x - delta,x + delta];            
            
        end % end of etaLogSigmaCI.
    
        % CIs on Canonical parameter vector.
        function CI = hetaSigmaCI(slme,alpha)
%hetaSigmaCI - Confidence intervals on the Canonical parameters.
%   CI = hetaSigmaCI(slme,alpha) takes a StandardLinearLikeMixedModel
%   object slme and a confidence level between 0 and 1 and computes the
%   (1-alpha)*100 CIs for elements of the vector [heta;sigma] where heta is
%   the Canonical parameter vector. CI(:,1) contains the lower limits of
%   the CIs and CI(:,2) contains the upper limits of the CIs. The last row
%   of CI contains the lower and upper limit of the CI for sigma.
        
            % (1) Ensure that alpha is sensible.
            assert( isnumeric(alpha) & isreal(alpha) & isscalar(alpha) );
            assert( alpha >= 0 & alpha <= 1 );

            % (2) Get CIs on [eta;log(sigma)].
            CI0 = etaLogSigmaCI(slme,alpha);
            
            % (3) Get the Psi object from slme.
            Psi = slme.Psi;
            
            % (3) CI0(:,1) contains lower bounds for a (1-alpha) CI on
            % [eta;log(sigma)]. Transform CI0(:,1) to a lower bound on
            % heta and sigma.
            eta_lb = CI0(1:end-1,1);
            logsigma_lb = CI0(end,1);
            sigma_lb = exp(logsigma_lb);
            Psi = setNaturalParameters(Psi,eta_lb);
            Psi = setSigma(Psi,sigma_lb);
            heta_lb = getCanonicalParameters(Psi);
            
            % (4) Transform CI0(:,2) to an upper bound on heta and sigma.
            eta_ub = CI0(1:end-1,2);
            logsigma_ub = CI0(end,2);
            sigma_ub = exp(logsigma_ub);
            Psi = setNaturalParameters(Psi,eta_ub);
            Psi = setSigma(Psi,sigma_ub);
            heta_ub = getCanonicalParameters(Psi);
            
            % (5) Assemble CI.
            CI(:,1) = [heta_lb;sigma_lb];
            CI(:,2) = [heta_ub;sigma_ub];                        
            
        end % end of hetaSigmaCI.                
    end        
    
    % Public utility methods.
    methods (Access=public)           
        function bsim = randomb(slme,S,NumReps)
%randomb - Generate random effects vectors using a fitted model.           
%   bsim = randomb(slme,S,NumReps) takes a fitted
%   StandardLinearLikeMixedModel slme and generates random effects vectors
%   in bsim from slme using the RandStream object S. If slme.Psi has K
%   blocks then NumReps is a length K vector. The vector bsim is a column
%   vector:
%
%   bsim = [b(1)_1;b(1)_2;...;b(1)_NumReps(1);...
%           b(2)_1;b(2)_2;...;b(2)_NumReps(2);...
%           b(K)_1;b(K)_2;...;b(K)_NumReps(K)]
%   where
%       b(i)_j = jth vector from ith matrix slme.Psi.Matrices{i}.
%
%   For block i, we generate NumReps(i) independent random vectors, each of
%   size Psi.SizeVec(i) from slme.Psi.Matrices{i} and concatenate them
%   vertically. We do this for all the blocks and return a concatenated
%   vector of size sum_{i = 1 to K} NumReps(i) * Psi.SizeVec(i).

            % (1) Extract Psi from slme. 
            Psi = slme.Psi;
                        
            % (2) Also extract SizeVec and NumBlocks from Psi.
            SizeVec = Psi.SizeVec;
            NumBlocks = Psi.NumBlocks;
            
            % (3) NumReps is a numeric, real, vector of length NumBlocks.
            assert( isnumeric(NumReps) & isreal(NumReps) ...
                & isvector(NumReps) & length(NumReps) == NumBlocks );
            
            % (4) Make NumReps into a column vector.
            if size(NumReps,1) == 1
                NumReps = NumReps';
            end
            
            % (5) Psi is a BlockedCovariance object. It has NumBlocks
            % smaller covariance matrices given in the Psi.Matrices cell
            % array. Psi.Matrices{i} has size SizeVec(i).
            
            % (6) The combined random effects vector bsim will be of size 
            % sum_{i=1 to NumBlocks} SizeVec(i)*NumReps(i).            
            bsim = zeros(sum(SizeVec .* NumReps),1);
                        
            % (7) Loop over NumBlocks. For block i, generate NumReps(i)
            % vectors from Psi.Matrices{i}, the length of each such vector
            % being SizeVec(i) times 1.
            if ~isempty(bsim)
                offset = 0;
                for r = 1:NumBlocks       
                    % (8) Get the covariance matrix PSIr for block r.
                    Lr = getLowerTriangularCholeskyFactor(Psi.Matrices{r});
                    sigmar = getSigma(Psi);
                    PSIr = (sigmar^2)*(Lr*Lr');           

                    % (9) Get the row vector for the mean of block r.
                    mur = zeros(1,SizeVec(r));

                    % (10) How many vectors to be drawn for the current block?
                    Nr = NumReps(r);

                    % (11) Generate random vectors in columns of temp.
                    if isempty(S)
                       % Use default mvnrnd to get a matrix with Nr columns.
                       temp = mvnrnd(mur,PSIr,Nr)';
                       % The following call should also give the same result:
                       % temp = slme.mymvnrnd([],mur,PSIr,Nr)';                                       
                    else
                       % Use mymvnrnd that accepts RandStream S. Get a matrix
                       % with Nr columns.
                       temp = slme.mymvnrnd(S,mur,PSIr,Nr)';
                    end

                    % (12) Concatenate temp to get a SizeVec(r)*NumReps(r) by 1
                    % vector and store in bsim starting from offset.
                    bsim(offset+1 : offset+SizeVec(r)*NumReps(r)) = temp(:);

                    % (13) Update offset for the next block.
                    offset = offset + SizeVec(r)*NumReps(r);                
                end
            end

        end % end of randomb.
        
        function slme = turnOffOptimizerDisplay(slme)
%turnOffOptimizerDisplay - Sets the display option to effective 'off'.            
%   slme = turnOffOptimizerDisplay(slme) takes a
%   StandardLinearLikeMixedModel object slme, and depending on the value of
%   slme.Optimizer, modifies the slme.OptimizerOptions property such that
%   subsequent calls to refit do not display the convergence information.
%   The modified object is returned.

            switch lower(slme.Optimizer)
                case 'quasinewton'
                    slme.OptimizerOptions.Display = 'off';
                case 'fminunc'
                    slme.OptimizerOptions.Display = 'off';
                case 'fminsearch'
                    slme.OptimizerOptions.Display = 'off';
                otherwise
                    error(message('stats:classreg:regr:lmeutils:StandardGeneralizedLinearMixedModel:BadOptimizerName'));
            end

        end % turnOffOptimizerDisplay.                
    end
    
    % Protected utility methods.
    methods (Access=protected)        
        function [thetaHat,cause] = doMinimization(slme,fun,theta0)
%doMinimization - Minimizes fun starting at theta0.
%   [thetaHat,cause] = doMinimization(slme,fun,theta0) attempts to minimize
%   fun starting at theta0 using the optimizer and optimizer options stored
%   in slme. This is a helper method for solveForThetaHat. The output
%   thetaHat is a best guess for the solution and cause is an integer that 
%   indicates the qualitative reason for termination of the optimizer. 
%   Here's what the integer codes in cause mean:
%
%   cause = 0      means that the problem was successfully solved. In most 
%                  cases this means that magnitude of the gradient is small 
%                  enough. The interpretation is "Local minimum found"
%                  (subject to Hessian checks).
%
%   cause = 1      means that problem may have been successfully solved. In
%                  most cases this means that the step size most recently 
%                  taken was small or the changes in objective function was
%                  small. The interpretation is "Local minimum possible".
%
%   cause = 2      means unable to converge to a solution. This may mean 
%                  that iteration/function evaluation limit was reached or 
%                  unable to converge to the specified tolerances.            

            switch lower(slme.Optimizer)                
                case 'quasinewton'                
                    % (1) Call the optimizer. The sequence of outputs is
                    % [theta,funtheta,gradfuntheta,cause].
                    [thetaHat,~,~,cause] = ...
                        fminqn(fun,theta0,'Options',slme.OptimizerOptions);
                    
                    % (2) For fminqn, cause already has the right meaning.
                    % So, no need for any translation.
                    
                case 'fminunc'                    
                    
                    % (1) Turn optim:fminunc:SwitchingMethod warning off.
                    % Since our default fminunc optimization algorithm is
                    % 'quasi-newton', we do not need to disable the 
                    % switching method warning. However, if we decide to
                    % make fminunc's 'trust-region' algorithm the default
                    % then we would need the next 3 lines.
                    %warnState = warning('query','all');
                    %warning('off','optim:fminunc:SwitchingMethod');
                    %cleanupObj = onCleanup(@() warning(warnState));
                    
                    % (2) Call the optimizer. The sequence of outputs is
                    % [X,FVAL,EXITFLAG].
                    [thetaHat,~,exitflag] = ...
                        fminunc(fun,theta0,slme.OptimizerOptions);
                    
                    % (3) Translate exitflag into the cause interpretation.
                    switch exitflag
                        case 1 
                            % Magnitude of gradient is small enough.
                            cause = 0;
                        case {2,3,5}
                            % Step of function tolerance reached.
                            cause = 1;
                        otherwise
                            % Unable to converge.
                            cause = 2;
                    end
                    
                case 'fminsearch'
                    % (1) Call the optimizer. The sequence of outputs is 
                    % [X,FVAL,EXITFLAG].
                    [thetaHat,~,exitflag] = ...
                        fminsearch(fun,theta0,slme.OptimizerOptions);
                    
                    % (2) Translate exitflag into the cause interpretation.
                    switch exitflag
                        case 1
                            % Maximum coordinate difference between current
                            % best point and other points in simplex is
                            % less than or equal to TolX, and corresponding
                            % difference in function values is less than or
                            % equal to TolFun.
                            cause = 0;
                        otherwise
                            % Unable to converge.
                            cause = 2;                            
                    end
                    
                otherwise
                    % <entry key="BadOptimizerName">''Optimizer'' must be either ''quasinewton'' or ''fminunc''.</entry>
                    error(message('stats:classreg:regr:lmeutils:StandardLinearMixedModel:BadOptimizerName'));
            end
            
        end % end of doMinimization.        
    end
    
    % Static utility methods.
    methods (Static)              
        function C = covarianceOnNaturalScale(H) 
%covarianceOnNaturalScale - Approximate the covariance of natural parameters.
%   C = covarianceOnNaturalScale(H) takes a symmetric matrix H where H is
%   either the Hessian of the negative profiled log likelihood (beta
%   profiled out) or the Hessian of the negative restricted log likelihood
%   as a function of [eta;log(sigma)] evaluated at [etaHat;log(sigmaHat)].
%   The covariance C is given by inv(H) except that here we account for
%   cases when H is not positive definite.

            % (1) Tentatively, set C as the pseudo inverse of H.
            C = pinv(H);

            % (2) Which elements of the Natural parameter vector are
            % estimable considering H as equivalent to the Normal matrix
            % X'*X in linear regression.
            sizeH = size(H,1);
            F = internal.stats.isEstimable(eye(sizeH),'NormalMatrix',H,'TolSVD',sqrt(eps(class(H))));

            % (3) If C(i,i) is < 0, F(i) must be set to false since the ith 
            % element of Natural parameter vector cannot have a negative
            % variance.
            diagC = diag(C);
            F(diagC < 0) = false;
            
            % (4) Set elements of non-estimable rows and columns as NaN.
            C(~F,:) = NaN;
            C(:,~F) = NaN;            
            
        end % end of covarianceOnNaturalScale.   
        
        function [ok,p,H,C] = fullrankH(H,C)
%FULLRANKH Make sure hypothesis matrix has full rank.
% [ok,p,H,C]=fullrankH(H,C) accepts a M by q matrix H and a M by 1 vector C
% used for testing the null hypothesis H*beta = C and returns a full row
% rank matrix H and a vector C that can be used for testing H*beta = C. ok
% is true if the input H and C are consistent otherwise it is false. p is 
% the rank of the input H.

% NOTE: The matrix H may not be of full row rank but H and C together must
% be sensible. For example, if beta is a 3 by 1 vector then H = [1 0 0] and
% C = [0] express the null hypothesis beta(1) = 0. The same hypothesis can
% be expressed by choosing H = [1 0 0;1 0 0] and C = [0 0] in which case
% the null hypothesis becomes [beta(1);beta(1)] = [0;0]. The second
% equation is redundant but the set of equations is still consistent. On
% the other hand, H = [1 0 0;1 0 0] and C = [1;0] express the null
% hypothesis [beta(1);beta(1)] = [1;0]. Clearly, these equations are not
% consistent.

            % (1) Find the rank of H and a basis for the row space.
            [~,r,e] = qr(H',0);
            tol = eps(norm(H))^(3/4);
            p = sum(abs(diag(r)) > max(size(H))*tol*abs(r(1,1)));

            % (2) Get a list of linearly dependent contrasts.
            E = e(1:max(1,p));

            % (3) Find coefficients that satisfy these contrasts.
            b = H(E,:)\C(E);

            % (4) Make sure they satisfy the entire set of contrasts.
            ok = all(abs(H*b-C) < tol);

            % (5) Now it is sufficient to use a full-rank subset.
            H = H(E,:);
            C = C(E);
    
        end % end of fullrankH.        
    
        function g = getGradient(fun,theta,step)
%getGradient - Get the numerical gradient by finite differences.
%   g = getGradient(fun,theta,step) gets the numerical gradient using
%   central differences of the function fun evaluated at p by 1 vector
%   theta. The output g will be a column vector of size p by 1. fun is a
%   function handle to a function that accepts a vector theta and returns a
%   scalar.

            % Set step size.
            if (nargin == 2)
                step = eps^(1/3);
            end

            % Initialize output.
            p = length(theta);
            g = zeros(p,1);   
            for i = 1:p   
                % Use central differences.
                theta1 = theta;
                theta1(i) = theta1(i) - step;

                theta2 = theta;
                theta2(i) = theta2(i) + step;

                g(i) = (fun(theta2) - fun(theta1))/2/step;        
            end
    
        end % end of getGradient.
        
        function H = getHessian2(fun,theta,wantRegularized)
%getHessian2 - Get the numerical Hessian by finite differences.
%   H = getHessian2(fun,theta,wantRegularized) gets the numerical Hessian by
%   central differences of the function fun evaluated at p by 1 vector
%   theta. The output H will be a p by p matrix. fun is a function handle
%   to a function that accepts a vector theta and returns a scalar. If
%   wantRegularized is true then we return H + deltaH where deltaH is given
%   by deltaH = delta*eye(size(H)), delta = abs(min(eig(H))) + sqrt(eps).
%   The idea is that H may have negative eigenvalues but H + deltaH is
%   positive definite. The default value of wantRegularized is false.

            % Set step size.
            step = eps^(1/4);

            % Initialize output.
            p = length(theta);    
            H = zeros(p,p);  %#ok<*PROP>
            for i = 1:p
                % Use central differences.
                theta1 = theta;
                theta1(i) = theta1(i) - step;

                theta2 = theta;
                theta2(i) = theta2(i) + step;

                H(:,i) = ( classreg.regr.lmeutils.StandardLinearLikeMixedModel.getGradient(fun, theta2, step) ...
                    -  classreg.regr.lmeutils.StandardLinearLikeMixedModel.getGradient(fun, theta1, step) )/2/step;       
            end

            if ( nargin == 3 && wantRegularized )        
                lambda = eig(H);
                if ~all(lambda > 0)
                    % There is a -ve or zero eigenvalue.
                    delta = abs(min(lambda)) + sqrt(eps);
                    % Add a multiple of identity to H.
                    H = H + delta*eye(size(H));
                end
            end
    
        end % end of getHessian2.
         
        function H = getHessian(fun, theta, wantRegularized)
% H = getHessian(fun, theta, wantRegularized) gets the numerical Hessian of
% the function fun evaluated at p-by-1 vector theta. The output H will be a
% p-by-p matrix. fun is a function handle to a function that accepts a
% vector theta and returns a scalar.
%
% Set the third argument to true to regularize the Hessian by adding a
% multiple of identity.

            % 1. Set step size.
            step = eps^(1/4);

            % 2. Initialize variables and evaluate fun(theta) once.
            p = length(theta);
            H = zeros(p,p);
            funtheta = fun(theta);
            denom = 4*(step^2);

            % 3. Numerical derivatives using central differences.
            for i = 1:p
                for j = i:p

                    if ( j == i )
                        % 4. Diagonal element.
                        theta2 = theta;
                        theta2(i) = theta2(i) + 2*step;

                        theta1 = theta;
                        theta1(i) = theta1(i) - 2*step;

                        H(i,j) = (fun(theta2) + fun(theta1) - 2*funtheta)/denom;
                    else
                        % 5. Non-diagonal element.
                        theta4 = theta;
                        theta4(i) = theta4(i) + step;
                        theta4(j) = theta4(j) + step;

                        theta3 = theta;
                        theta3(i) = theta3(i) - step;
                        theta3(j) = theta3(j) + step;

                        theta2 = theta;
                        theta2(i) = theta2(i) + step;
                        theta2(j) = theta2(j) - step;

                        theta1 = theta;
                        theta1(i) = theta1(i) - step;
                        theta1(j) = theta1(j) - step;

                        H(i,j) = (fun(theta4) + fun(theta1) - fun(theta3) - fun(theta2))/denom;
                    end
                end
            end
    
            % 6. Fill in lower triangular part.
            H = triu(H,1)' + H;    

            % 7. Regularize if required.
            if ( nargin == 3 && wantRegularized )        
                lambda = eig(H);
                if ~all(lambda > 0)
                    % There is a -ve or zero eigenvalue.
                    delta = abs(min(lambda)) + sqrt(eps);
                    % Add a multiple of identity to H.
                    H = H + delta*eye(size(H));
                end
            end
            
        end % end of getHessian.        
        
        function R = mymvnrnd(S, MU, SIGMA, N)
% This helper function works exactly like mvnrnd, except that it uses a
% specified RandStream object to generate samples.
%
% INPUTS:
%
%     S = A RandStream object from which to draw samples. This can be empty
%         in which case, the default RandStream will be used.
%    MU = 1 by p mean vector of the multivariate Normal distribution.
% SIGMA = p by p co-variance matrix of the multivariate Normal distribution.
%     N = number of 1 by p vectors to sample from the multivariate Normal distribution.
%
% OUTPUTS:
%     R = N by p matrix where each row is a draw from a multivariate Normal
%         distribution with mean vector MU and co-variance matrix SIGMA.
% 
%
% Suppose x is a p by 1 vector with x ~ N(mu, sigma) where mu is a p by 1
% mean vector and sigma is a p by p co-variance matrix. The function cholcov 
% computes a factor T such that:
% 
% T = cholcov(sigma) and 
% sigma = T'*T.
%
% Now,
%        x ~ N(mu, sigma) 
%
%       x - mu ~ N(0, sigma)
%
%     inv(T')*(x-mu) ~ N(0, I_p)
%
% Suppose y ~ N(0, I_p) is a p by 1 vector (using randn for example) then
%
% inv(T')*(x-mu) = y
%
% Thus,  x = mu + T'*y
%
%        x' = mu' + y' * T
%
% If we want N vectors x_1,x_2,..,x_N from N(mu, sigma) then:
%
%       [x_1'       [mu'     [y_1'
%        x_2'        mu'      y_2'
%        x_3'    =   mu'  +   y_3'   * T
%         .           .        .
%         .           .        .
%        x_N']       mu']     y_N']
%
% The N by p matrix formed from y_1,y_2...,y_N can be generated using
% randn(N,p). Also note that the input to our function MU = mu'. Hence we
% can generate the N by p matrix made up of x_1,...,x_N (which is the same
% as the output of our function R) as follows;
%
%       [x_1'       
%        x_2'        
%        x_3'    =  R = N by p matrix = repmat(MU,N,1) + randn(N,p) * T; 
%         .           
%         .           
%        x_N']      
%

            % Get p and make sure MU and SIGMA are of the right size.
            p = size(MU,2);
            assert(p == size(SIGMA,1));
            assert(p == size(SIGMA,2));

            % Compute the factor T.
            T = cholcov(SIGMA);

            % Compute the output matrix R.
            if isempty(S)
                % Use default RandStream.
                R = bsxfun(@plus, MU, randn(N,p)*T);
            else
                % Use user supplied RandStream.
                assert( isa(S,'RandStream') );
                R = bsxfun(@plus, MU, randn(S,N,p)*T);
            end
    
        end % end of mymvnrnd.
        
        function tf = hasNaNInf(X)
%hasNaNInf - Check if input has NaNs or Infs.
%   tf = hasNaNInf(X) takes a matrix/vector X and if X has NaNs or Infs 
%   anywhere in it then tf is true otherwise tf is false.

            tf = any(any(isnan(X))) || any(any(isinf(X)));

        end % end of hasNaNInf. 
        
        function assertThat(condition,msgID,varargin)
%assertThat(condition,msgID,varargin) takes a variable condition that is
%   either true or false, a message catalog ID msgID and optional arguments
%   required to create a message object from msgID. If condition is false,
%   an error message represented by msgID is thrown.

            if ~condition        
                % (1) Create a message object from msgID and varargin.
                try
                    msg = message(msgID,varargin{:});
                catch
                    % <entry key="BadMsgID">Invalid message ID: {0}.</entry>      
                    error(message('stats:LinearMixedModel:BadMsgID',msgID));
                end        
                % (2) Create and throw an MException.
                ME = MException( msg.Identifier, getString(msg) );
                throwAsCaller(ME);        
            end

        end % end of assertThat.
        
        function x = logAbsDetTriangular(R)
%logAbsDetTriangular - Log absolute determinant of triangular matrix.
%   x = logAbsDetTriangular(R) takes an upper triangular or lower
%   triangular matrix R and computes the output x = log( abs(det(R)) ). R
%   is assumed to be triangular, we do not explicitly check this in this
%   method. If R is not triangular, the result will be incorrect.

%   Note: If R is of size p by p then det(R) = r_11 * r_22 * ... * r_pp and
%   so abs(det(R)) = abs(r_11) * abs(r_22) * ... * abs(r_pp). Therefore,
%   log( abs(det(R)) ) 
%       = log( abs(r_11) * abs(r_22) * ... * abs(r_pp) )
%       = log( abs(r_11) ) + log( abs(r_22) ) + ... + log( abs(r_pp) )
%       = sum_{i=1}^p log( abs(r_ii) )
%       = sum(log(abs( diag(R) )));

            x = sum(log(abs( diag(R) )));
    
        end % end of logAbsDetTriangular.  
                
        function checkPositiveDefinite(H,msgID,varargin)
%checkPositiveDefinite - Checks if the input matrix is positive definite.
%   checkPositiveDefinite(H,msgID,varargin) takes a square matrix H that is
%   expected to be symmetric positive definite but may have negative
%   eigenvalues due to numerical roundoff. If a Cholesky factorization of H
%   fails, a warning is displayed to the user with the content given by the
%   message ID msgID. varargin can contain additional arguments required to
%   construct a message with ID msgID.
         
            try 
                chol(H);
            catch ME %#ok<NASGU> 
                % (1) Construct a message object and extract its string content.
                msg = message(msgID,varargin{:});
                msgStr = getString(msg);            

                % (2) Does H have any NaN's or Inf's?
                hasNaN = any(isnan(H(:)));
                hasInf = any(isinf(H(:)));           

                % (3) Compute the smallest eigenvalue of H if it has no NaN or
                % Inf values. Otherwise, say that H has Inf/NaN values.
                if ~hasNaN && ~hasInf
                    mineigH = min(eig(H));
                    % <entry key="Message_MinEig">Minimum eigenvalue is {0}.</entry>
                    eigmsgStr = getString(message('stats:classreg:regr:lmeutils:StandardLinearLikeMixedModel:Message_MinEig',num2str(mineigH)));
                    msgStr = [msgStr,' ',eigmsgStr];
                else
                    % <entry key="Message_NaNInfInHessian">Hessian has NaNs or Infs.</entry>
                    naninfmsgStr = getString(message('stats:classreg:regr:lmeutils:StandardLinearLikeMixedModel:Message_NaNInfInHessian'));
                    msgStr = [msgStr,' ',naninfmsgStr];
                end
                % (4) Throw a warning with the message ID = msg.Identifier.
                warning(msg.Identifier,msgStr);
            end
        
        end % end checkPositiveDefinite.  
        
        function mu = constrainVector(mu,a,b)
%mu = constrainVector(mu,a,b) takes a vector mu and scalars a and b with
%   a < b and constrains mu between a and b. Values of mu < a are set equal
%   to a and values of mu > b are set equal to b.

            assert(a < b);
            mu(mu < a) = a;
            mu(mu > b) = b;
    
        end % end of constrainVector.    
        
        function p = solveUsingQR(J,r)
%p = solveUsingQR(J,r) takes a q-by-q matrix J, a q-by-1 vector r and
%   computes the solution p that minimizes the 2-norm of (J*p + r).
%
%%
% We want to solve $0.5 ||J p + r||_2^2 \rightarrow min$ or a solution to
% the equation:
%
% $J^T J \, p = -J^T r$
%
% The matrix $J$ is decompsed via a $QR$ factorization as follows:
%
% $JE = QR$
%
% where $E$ is a permutation matrix satisfying 
% 
% $EE^T = I_q$. 
% 
% Let $diagR$ be the diagonal of the $R$ matrix. We find an index vector
% $idx$ such that $diagR(idx)$ are non-zero (numerically speaking). For
% this purpose, we can use a test like this:
%
% $tol = max(abs(diagR))*sqrt(eps(class(diagR))))$
%
% $idx = find(abs(diagR) > tol)$
%
% If $idx$ is empty, we return a vector of all zeros as the solution $p$.
% Otherwise, we partition matrices $Q$ and $R$ using $idx$ as follows:
%
% $Q_1 = Q(:,idx)$ and
% $R_{11} = R(idx,idx)$
%
% Let $p_{E} = zeros(q,1)$. Suppose $length(idx) = m$ then we can partition
% $p_{E} = [p_{E1};p_{E2}]$ where $p_{E1}$ is $m \times 1$ and $p_{E2}$ is
% $(q-m) \times 1$. 
%
% $p_{E1}$ is given by: $p_{E1} = -R_{11}^{-1} Q_1^T r$
%
% and $p_{E2}$ is given by: $p_{E2} = zeros(q-m,1)$
%
% From $p_{E}$, the solution $p$ is obtained as follows:
%
% $p = zeros(q,1)$
%
% $E^T p = p_{E}$
%
% If the permutation information is returned as a vector $e$ instead of a
% matrix $E$ then:
%
% $p(e) = p_{E}$.

            % 1. Validate inputs.
            assert(isnumeric(J) & isreal(J) & ismatrix(J));
            assert(isnumeric(r) & isreal(r) & iscolumn(r));
            q = size(J,1);
            assert(size(r,1) == q);

            % 2. QR factorization of J.
            [Q,R,e] = qr(J,'vector');

            % 3. Get index of non-zero elements on diagonal of R.
            diagR = diag(R);
              tol = max(abs(diagR))*sqrt(eps(class(diagR)));
              idx = find(abs(diagR) > tol);

            % 4. Decide what to do based on idx.
            if isempty(idx)
                % 4.1 If idx is empty we are done.
                p = zeros(q,1);
            else        
                % 4.2. Get Q1 and R11.
                 Q1 = Q(:,idx);
                R11 = R(idx,idx);

                % 4.3 Get pE.
                pE1 = -(R11 \ (Q1'*r));
                  m = length(idx);
                pE2 = zeros(q-m,1);
                 pE = [pE1;pE2];

                % 4.4 Get p.
                   p = zeros(q,1);
                p(e) = pE;
            end
        
        end % end of solveUsingQR.        
    end
        
end

