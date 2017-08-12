function [theta,funtheta,gradfuntheta,cause] = fminqn(fun,theta0,varargin)
%

%fminqn Utility to solve unconstrained minimization problems.
%[theta,funtheta,gradfuntheta,cause] = fminqn(fun,theta0,varargin)
% solves the minimization problem:
%
%             min w.r.t theta: fun(theta)
%
% starting from the initial point theta0. A trust region based quasi-Newton
% method using the symmetric rank-1 (SR1) approximation to the Hessian is
% used. Trust region sub-problems are solved using Steihaug's modified
% conjugate gradient (CG) method. If negative curvature is encountered
% during Steihaug's CG, then the trust region problem is solved exactly.
% Gradient of fun is computed using finite differences. The function fun is
% assumed to be smooth.
%
%   POSITIONAL PARAMETERS:
%
%           fun   A function handle to the function to be minimized that 
%                 can be called like fun(theta0). 
%
%        theta0   Initial point to begin iterations (column vector).
%
%
%   OPTIONAL NAME/VALUE PAIRS:
%
%       Options   A structure containing optimization options with the 
%                 following fields:
% 
%          'TolFun'  -  Relative tolerance on the gradient of the objective
%                       function. Default is 1e-6. Suppose f and gradf are
%                       the values of fun and gradient of fun evaluated at
%                       the current iterate. Also let gradf0 be the
%                       gradient of fun evaluated at the initial point
%                       theta0. Iterations stop if:
%
%   max(abs(gradf)) <= max( 1, min( abs(f), max(abs(gradf0)) ) ) * TolFun 
%                    
%             'TolX'  - Absolute tolerance on the step size. Default is 
%                       1e-12. If s is the current step then iterations 
%                       stop if:
%
%                           norm(s) <= TolX
%  
%          'MaxIter'  - Maximum number of iterations allowed. Default is
%                       10000.
%  
%          'Display'  - Level of display:  'off', 'iter', or 'final'.
%                       Default is 'off'. Both 'iter' and 'final' have the
%                       same effect.
%
%   OUTPUTS:
%
%         theta    Estimated solution to the problem.
%
%      funtheta    Function value at the optimal theta.
%
%  gradfuntheta    Gradient of the function at optimal theta.
%
%         cause    An integer indicating the reason for termination.
%                  Possible values are as follows:
%                 
%                  cause = 0 means that the maximum absolute gradient of
%                  the function was less than TolFun in a relative error
%                  sense (see above for more details). The interpretation
%                  is: "Local minimum found".
%
%                  cause = 1 means that the step size most recently taken
%                  was less than TolX in an absolute error sense. The
%                  interpretation is "Local minimum possible".
%
%                  cause = 2 means 'iteration limit reached'.
%                   
% NOTES:
%
%   (1) What if we want to maximize instead of minimize? 
%   We normally minimize fun(theta). If you want to maximize fun(theta),
%   create a new function h(theta) = -fun(theta) and then minimize
%   h(theta).
%
%   (2) How to account for bound constraints on theta?
%   Suppose theta is a vector that satisfies a <= theta <= b where a and b
%   are the lower and upper bound vectors for elements of theta. We can
%   eliminate the bounds via the transformation:
%
%           theta = a + (b - a)./(1 + exp(-\phi))
%
%   Our objective function is now parameterized by phi and can be written
%   as:
%
%           h(phi) = fun(theta) = fun( a + (b - a)./(1 + exp(-\phi)) )
%
%   For any phi, we have a <= theta <= b. Hence, minimizing h(phi) w.r.t
%   phi is the same as minimizing fun(theta) while respecting the bound
%   constraints.
%
%   (3) It is the responsibility of the caller to ensure that fun does not
%   return NaN/Inf values.

%   Copyright 2012-2013 The MathWorks, Inc.

    %% Handle input args.        
    
        % (1) Check number of input arguments.
        narginchk(2,Inf);
    
        % (2) Extract optimization options.
        
            % (1) An 'Options' structure with default values for TolFun,
            % TolX, Display and MaxIter.
            dfltTolFun = 1e-6;
            dfltTolX = 1e-12;
            dfltDisplay = 'off';
            dfltMaxIter = 10000;            
            dfltoptions = ...
                statset('TolFun',dfltTolFun,...
                'TolX',dfltTolX,...
                'Display',dfltDisplay,...
                'MaxIter',dfltMaxIter);            
            % (2) Process optional name/value pairs.                        
                names = {'Options'};
                dflts = {dfltoptions};
                options = ...
                    internal.stats.parseArgs(names,dflts,varargin{:});                
            % (3) Validate options.                   
                if ~isstruct(options)
                    % <entry key="BadOptions">The ''Options'' input must be a structure.</entry>
                    error(message('stats:classreg:regr:lmeutils:fminqn:BadOptions'));
                end                
            % (4) Combine dfltoptions and user specified options to create
            % a structure options with values for TolFun, TolX, Display
            % and MaxIter.
                options = statset(dfltoptions,options);        
            % (5) Extract gradTol, stepTol and maxit from the options
            % structure.
                gradTol = options.TolFun;
                stepTol = options.TolX;
                  maxit = options.MaxIter;                
            % (6) If options.Display is 'off' then set verbose to false,
            % otherwise set it to true.
                if strcmpi(options.Display,'off')
                    verbose = false;
                else
                    verbose = true;
                end

        % (3) Validate fun and theta0.
        
            % (1) fun must be a function handle.
            fun = validateFun(fun);            
            
            % (2) theta0 must be a numeric, real vector. If theta0 is a row
            % vector, we will convert it into a column vector. theta0 is
            % not allowed to contain NaN/Inf values.
            theta0 = validateTheta0(theta0);            
            
    %% Edge case when theta0 is empty.
    
        if ( numel(theta0) == 0 )            
                   theta = theta0;
                funtheta = [];
            gradfuntheta = [];
                   cause = 0;
            return;
        end
        
    %% Initialize control variables for trust region iterations.  
    
        % (1) Make function to compute gradient. Instead of doing 
        % getGradient(fun,x), we can now do gradfun(x).
        gradfun = MakeGradient(fun);
    
        % (2) Initial value of solution, initial function and gradient. 
        % We will also save the infinity norm of g at every iteration. We 
        % will always keep x, f, g and infnormg syncronized. We guard
        % against +Inf or -Inf values in g.
                   x = theta0;
                   f = fun(x);
                   errorIfNotScalar(f);
                   g = gradfun(x);        
                   g = replaceInf(g,realmax);
            infnormg = max(abs(g));
                    
        % (3) Edge case when f is -Inf.
        if ( f == -Inf )
            % Can't do any better than this.
                   theta = x;
                funtheta = f;
            gradfuntheta = g;
                   cause = 0;
            return;
        end
            
        % (4) Initial positive definite Hessian approximation. For problems
        % with more than roughly 20 variables, we don't want to wast time
        % in computing the full numerical Hessian and we settle for a
        % diagonal Hessian as a reasonable approximation for B.
        if ( length(x) > 21 )
            B = getDiagonalHessian(fun,x);
            Imult = max(100,max(abs(diag(B))));
            B = Imult*eye(length(x));
        else
            B = getHessian(fun,x);
        end

        % (5) Guard against Inf or -Inf values in B.    
        B = replaceInf(B,realmax);
        
        % (6) Save initial infinity norm of gradient (i.e., at theta0) for 
        % relative gradient convergence test.
        infnormg0 = infnormg;             
        
        % (7) Initial trust region radius. Other choices are possible. For 
        % non-diagonal B, we may do: Delta = max(1,max(abs(pinv(B)*g))).
        
            % (1) Overall bound on step lengths.
            DeltaMax = 1e9;
            
            % (2) Initialize Delta based on whether B is positive definite 
            % or not.            
            [p0,~,isposdefB] = solveTrustRegionProblemExact(g,B,DeltaMax);            
            if isposdefB
               Delta = max(sqrt(length(x)), norm(p0));
            else
               Delta = sqrt(length(x));                                         
            end 
               
        % (8) Sufficient decrease and skipping condition for SR1 method.
        % It is unlikely that these parameters will need to be changed.
            eta = 5e-4;
              r = 1e-8;
      
        % (9) Convergence flag. The while loop that follows will stop if
        % found = true. Iterations can stop either when optimal solution is
        % found or when step size becomes smaller than the specified
        % tolerance.
        found = false;
            
        % (10) Iteration counter starting from 0. Do not change this.
        iter = 0;                
        
        % (11) Number of steps accepted so far.
        numAccepted = 0;
        
    %% Solve a sequence of trust region subproblems and change size of trust region.    
        while( found == false )        

            % (1) Solve trust region problem using current values of g, B
            % and Delta. This will give us a step size s.

                % (1) Convergence tolerance for Trust region subproblem.
                epsk = min(0.5,sqrt(infnormg))*infnormg;
                [s,reasonCGTerm] ...
                    = solveTrustRegionProblem(g, B, Delta, epsk);
                if any(strcmpi(reasonCGTerm,{'NEG CURV','NaNORInf'}))   
                     % (2) Exact solution to the trust region problem.
                     B = replaceInf((B+B')/2,realmax);
                     s = solveTrustRegionProblemExact(g, B, Delta);
                     reasonCGTerm = 'EXACT';
                end

            % (2) Tentative new value of x.
            xs = (x+s);
                
            % (3) Compute y, the difference in gradients between x and
            % (x+s). The gradient at x is already known and equal to g.
                gs = gradfun(xs);
                 y = gs - g; 

            % (4) Compute actual and predicted reduction in fun.
                  fs = fun(xs);
                ared = f - fs;
                pred = -(g'*s + 0.5* (s'*B*s));

            % (5) Ratio of actual reduction to predicted reduction.
            rho = ared/pred;
            
            % (6) Update x, f, g and infnormg if we achieve sufficient
            % reduction.
                stepTaken = false;
                if (rho > eta)
                 numAccepted = numAccepted + 1;
                   stepTaken = true;
                           x = xs;
                           f = fs; % also update the function at current x.
                           g = gs; % also update the gradient at current x.
                    infnormg = max(abs(g));
                end
            
            % (7) Get 2 norm of s.
            twonorms = norm(s);
            
            % (8) Display convergence info if requested.
                if (verbose == true)
                    displayConvergenceInfo(iter, f, infnormg, twonorms, reasonCGTerm, rho, Delta, stepTaken);
                end
            
            % (9) Update trust region radius.
                if (rho > 0.75)

                    if (twonorms > 0.8*Delta)
                        % Increase Delta but at most DeltaMax.
                        Delta = min(2*Delta,DeltaMax);
                    end

                elseif (rho < 0.1)

                    % Decrease Delta.
                    Delta = 0.5*Delta;

                end 

            % (10) Check skipping condition.
                y_Bs = (y - B*s);        
                y_Bs_times_s = s'*y_Bs;
                applyUpdate = abs(y_Bs_times_s) >= r * twonorms * norm(y_Bs);                            
                
            % (11) Update B.
                if applyUpdate
                   % (1) Tentative increment for B.
                   incr = ( (y_Bs*y_Bs')/y_Bs_times_s );
                   if ~any(isnan(incr(:))) && ~any(isinf(incr(:)))
                       % (2) No NaNs in incr. Increment B. 
                       B = B + incr;                        
                   else
                       % Warn the user about NaNs in increment.
                       % warning('MATLAB:NaNSR1Update','Increment to SR1 Hessian approximation had NaNs. Skipping update');
                   end
                end            
                
            % (12) Check convergence using relative error.
                tau = max(1, min( abs(f), infnormg0 ));   
                if (infnormg <= tau * gradTol)
                    found = true;
                    % Local minimum found.
                    cause = 0;
                elseif (twonorms <= stepTol)
                    found = true;
                    % Local minimum possible.
                    cause = 1;    
                elseif (f == -Inf)
                    % Can't do any better.
                    found = true;
                    % Local minimum found.
                    cause = 0;
                end            
            
            % (13) Check iteration counter.
                if ( iter > maxit )
                    found = true;
                    % Itetation limit reached.
                    cause = 2;
                end
                
            % (14) Display final convergence message.
                if (found == true && verbose == true)
                    displayFinalConvergenceMessage(infnormg, tau, gradTol, twonorms, stepTol, cause);  
                end
                
            % (15) Update iteration counter.
            iter = iter + 1;
                
        end % end of while.

        %% Set outputs.
            theta = x;
            if nargout > 1
                   funtheta = fun(x); 
               gradfuntheta = gradfun(x);
            end
    
end % end of fminqn.

%=== displayFinalConvergenceMessage
function displayFinalConvergenceMessage(infnormg, tau, gradTol, twonorms, stepTol, cause)
% Helper function to display final convergence message.
%      
%  infnormg = infinity norm of the gradient at current iterate.
%       tau = relative convergence factor for relative convergence test.
%   gradTol = relative tolerance for convergence test. 
%  twonorms = 2 norm of the proposed step size s.
%   stepTol = absolute tolerance on step length.
%     cause = reason for termination of the algorithm.
%
    % Explain why iterations stopped.
    fprintf('\n');
    infnormgStr = getString(message('stats:classreg:regr:lmeutils:fminqn:FinalConvergenceMessage_InfNormGrad'));
    fprintf(['         ',infnormgStr,' ','%6.3e\n'], infnormg);
    twonormsStr = getString(message('stats:classreg:regr:lmeutils:fminqn:FinalConvergenceMessage_TwoNormStep'));
   	fprintf(['              ',twonormsStr,' ','%6.3e, ','TolX   =',' ','%6.3e\n'], twonorms, stepTol);
    relinfnormgStr = getString(message('stats:classreg:regr:lmeutils:fminqn:FinalConvergenceMessage_RelInfNormGrad'));
    fprintf([relinfnormgStr,' ','%6.3e, ','TolFun =',' ','%6.3e\n'], infnormg/tau, gradTol);

    % Explain what the final solution means.
    if ( cause == 0 )
        % Local minimum found.      
        fprintf([getString(message('stats:classreg:regr:lmeutils:fminqn:Message_LocalMinFound')),'\n']);
        %fprintf('Local minimum found.\n');
    elseif ( cause == 1 )
        % Local minimum possible.
        fprintf([getString(message('stats:classreg:regr:lmeutils:fminqn:Message_LocalMinPossible')),'\n']);
        %fprintf('Local minimum possible.\n');
    elseif( cause == 2 )
        % Iteration limit reached.
        fprintf([getString(message('stats:classreg:regr:lmeutils:fminqn:Message_UnableToConverge')),'\n']);
        %fprintf('Unable to converge to a solution.\n');
    end

end % end of displayFinalConvergenceMessage.

%=== displayConvergenceInfo
function displayConvergenceInfo(iter, f, infnormg, twonorms, reasonCGTerm, rho, Delta, stepTaken)
% Helper function to display iteration wise convergence info.
%
%         iter = iteration number.
%            f = function value at current iterate.
%     infnormg = infinity norm of the gradient at current iterate.
%     twonorms = infinity norm of the proposed step size s.
% reasonCGTerm = reason why CG was terminated. One of 'CONV', 'NEG CURV' or 'BNDRY'.
%          rho = ratio of actual reduction to predicted reduction.
%        Delta = current size of trust region.
%    stepTaken = true if step was taken and otherwise false.
%
%  We will display convergence info like this:
%  ===================================================================================
%  ITER    FUN VALUE    NORM GRAD    NORM STEP    CG TERM    RHO    TRUST RAD   ACCEPT
%  ===================================================================================
%     0    4.395e-01    3.501e-01    2.611e-03    CONV       0.34   4.633e-03     YES
%     1    4.388e-01    6.802e-03    2.092e-03    CONV       0.127  4.633e-03     YES
%     2    4.388e-01    3.320e-03    3.389e-03    BNDRY     -0.34   4.633e-03     YES
%     3    4.388e-01    8.652e-04    1.383e-03    NEG CURV   0.83   2.316e-03     YES
%

    if ( rem(iter,20) == 0 )
        % Display header.
        fprintf('\n');
        fprintf('  ============================================================================================\n');
        fprintf('  ITER     FUN VALUE    NORM GRAD    NORM STEP    CG TERM        RHO        TRUST RAD   ACCEPT\n');
        fprintf('  ============================================================================================\n');
    end
    
    % Display iteration wise convergence info.
    if stepTaken == true
        stepTakenString = 'YES';
    else
        stepTakenString = ' NO';
    end
    fprintf('%6d    %+6.3e    %06.3e    %06.3e    %8s    %+6.3e    %06.3e     %3s\n', iter, f, infnormg, twonorms, reasonCGTerm, rho, Delta, stepTakenString);    
end % end of displayConvergenceInfo.

%=== solveTrustRegionProblem
function [p, reasonCGTerm] = solveTrustRegionProblem(g, B, Delta, epsk)
% Helper function to solve trust region sub-problems.
%
% INPUTS:
%
%     g = gradient of function fun at current iterate x.
%     B = quasi-Newton approximation to Hessian at current iterate x.
% Delta = current trust region radius.
%  epsk = convergence tolerance.
%
% We form a quadratic approximation to the function fun around the current
% iterate x and approximately solve the problem:
%
% min_{w.r.t p} m(p) = f + g'*p + 0.5*p'*B*p
%
% s.t. ||p||_2 <= Delta
%
% We use the CG-Steihaug algorithm to solve this problem.
%
% We do not need to know the value of fun at x (denoted by f above) to
% solve this problem and so f is not an input to this function. If g or B
% contains NaN or Inf we simply return a vector of NaNs.
%
% OUTPUTS:
%   
%                p = estimated solution vector.
%     reasonCGTerm = reason why CG iterations were terminated.
%
%     reasonCGTerm = 'CONV' means solution was found to desired tolerance.
%                  = 'BNDRY' means trust region boundary was reached.
%                  = 'NEG CURV' means negative curvature was detected.
%                  = 'NaNORInf' means NaN or Inf values were detected.
%

    % Protect against NaN and Inf.
    if ( any(isnan(g)) || any(isnan(B(:))) ...
            || any(isinf(g)) || any(isinf(B(:))) )
        p = NaN(length(g),1);
        reasonCGTerm = 'NaNORInf';
        return;
    end

    % Initialization.
    
        % Initial values of z, r and d.
            z = zeros(length(g),1);
            r = g;
            d = -r;

        % Why was CG terminated?
        reasonCGTerm = 'CONV';
            
        % Initial convergence check.
            if ( norm(r) < epsk )
                p = z;
                reasonCGTerm = 'CONV';
                return;
            end % end of if.

        % Convergence flag. We will stop iterating when found = true.
        found = false;
    
    % Begin steihaug's algorithm.    
        while( found == false )

           dBd = (d'*B*d);

           if ( dBd <= 0 )

              % Reason CG Terminated.
              reasonCGTerm = 'NEG CURV';
               
              % Find tau, such that p = z + tau*d minimizes m(p) and satisfies
              % ||p||_2 = Delta.
              %
              % If we solve || z + tau*d ||_2 = Delta, we get the following 2
              % values of tau:
              %
              % tau = [ -(z'*d) +- sqrt( (z'*d)^2 + (d'*d) * (Delta^2 - z'*z) ) ]/(d'*d);
              %
              % The quadratic equation that needs to be solved is:
              %
              %  a * tau^2 + b * tau + c = 0 where
              %
              %  a = d'*d;
              %
              %  b = 2*(z'*d);
              %
              %  c = (z'*z - Delta^2) 
              %
                  % Store various products for computing tau1 and tau2.
                      zd = z'*d;
                      dd = d'*d;
                      zz = z'*z;

                  % Solve quadratic equation and avoid cancellation error.    
                      [tau1, tau2] = solveQuadraticEquation( dd, 2*zd, (zz - Delta^2) );
                    %                   % Square root of Discriminant.
                    %                   sqrtDisc = sqrt( zd^2 + dd * (Delta^2 - zz) );
                    % 
                    %                   % Possible values of tau.
                    %                       tau1 = ( -zd + sqrtDisc )/dd;
                    %                       tau2 = ( -zd - sqrtDisc )/dd;

                  % Values of m(p) at tau1 and tau2.
                      p1 = z + tau1*d;
                      f1 = g'*p1 + 0.5*(p1'*B*p1);

                      p2 = z + tau2*d;
                      f2 = g'*p2 + 0.5*(p2'*B*p2);

                  % Choose p1 if f1 is smaller and choose p2 if f2 is smaller.  
                      if ( f1 <= f2 )
                          p = p1;
                      else
                          p = p2;
                      end % end of if.
                      
              return;

           end % end of if.

           alpha = (r'*r)/dBd;
            znew = z + (alpha*d);

           if ( norm(znew) >= Delta )

              % Reason CG Terminated.
              reasonCGTerm = 'BNDRY'; 
               
              % Find tau >= 0 such that p = z + tau*d satisfies ||p||_2 =
              % Delta. We proceed like before (see the case when dBd <= 0)
              %
              % The quadratic equation that needs to be solved is:
              %
              %  a * tau^2 + b * tau + c = 0 where
              %
              %  a = d'*d;
              %
              %  b = 2*(z'*d);
              %
              %  c = (z'*z - Delta^2) 
              %
              
              % Store various products for computing tau1 and tau2.
                  zd = z'*d;
                  dd = d'*d;
                  zz = z'*z;
                  
              % Solve quadratic equation and avoid cancellation error.
                  [tau1, tau2] = solveQuadraticEquation( dd, 2*zd, (zz - Delta^2) );                  
                %               % Square root of Discriminant.
                %               sqrtDisc = sqrt( zd^2 + dd * (Delta^2 - zz) );
                % 
                %               % Possible values of tau.
                %                   tau1 = ( -zd + sqrtDisc )/dd;
                %                   tau2 = ( -zd - sqrtDisc )/dd;

              % If tau1 is >=0 then set p = z + tau1*d and vice versa.
                  if ( tau1 >= 0 )
                      p = z + tau1*d;
                  else
                      p = z + tau2*d;
                  end % end of if.

              return;

           end % end of if.

           rnew = r + alpha*(B*d);

           if ( norm(rnew) < epsk )
              p = znew;
              reasonCGTerm = 'CONV';
              return;
           end % end of if.

           if any(isnan(rnew))
               p = NaN(length(g),1);
               reasonCGTerm = 'NaNORInf';
               return;
           end
           
           betanew = (rnew'*rnew)/(r'*r);       
              dnew = -rnew + betanew*d;       

           % Update variables for the next iteration.
           z = znew;
           r = rnew;
           d = dnew;
                    
        end % end of while.
    
end % end of solveTrustRegionProblem.

%=== solveQuadraticEquation
function [tau1,tau2] = solveQuadraticEquation(a,b,c)
% Find tau such that:
%
% a*tau^2 + b*tau + c = 0
%
% (1) Ensure that the discriminant D = b^2 - 4*a*c is >= 0. Otherwise throw a warning and
% return NaNs for tau1 and tau2.
%
% (2) If b > 0, then avoid cancellation error by defining:
%
%   tau1 = -1*(b + sqrt(D))/(2*a); and
%
%   tau2 = (-2*c)/(b + sqrt(D))
%
% (3) If b <= 0 then -b >= 0. Define tau1 and tau2 like this:
%
%   tau1 = (2*c)/(-b + sqrt(D));
%
%   tau2 = (-b + sqrt(D))/(2*a);
%

    % Compute discriminant.
    D = (b^2 - 4*a*c);

    % Make sure D is real and positive.
    assert( isreal(D) && D >= 0 );

    % Avoid cancellation errors.
    if (b > 0)
       
        b_plus_sqrtD = b + sqrt(D);
        
        tau1 = -1*b_plus_sqrtD/(2*a);
        
        tau2 = (-2*c)/b_plus_sqrtD;
        
    else
        
        minusb_plus_sqrtD = -b + sqrt(D);
        
        tau1 = (2*c)/minusb_plus_sqrtD;
        
        tau2 = minusb_plus_sqrtD/(2*a);
        
    end
    
end % end of solveQuadraticEquation.

%=== MakeGradient.
function g = MakeGradient(fun)

    g = @(Theta) getGradient(fun, Theta);

end % end of MakeGradient. 

% %=== MakeHessian. 
% function H = MakeHessian(fun)
% 
%     H = @(Theta) getHessian(fun, Theta);
% 
% end % end of MakeHessian.

%=== getGradient
function g = getGradient(fun, theta, step)
% Get numerical gradient of the function fun evaluated at p by 1 vector theta.
% The output g will be a column vector of size p by 1. fun is a function handle 
% to a function that accepts a vector theta and returns a scalar.

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
    
    g = replaceInf(g,realmax);
    
end % end of getGradient.

%=== getHessian
% function H = getHessian(fun, theta, wantRegularized)
% % Get numerical Hessian of the function fun evaluated at p by 1 vector theta.
% % The output H will be a p by p matrix. fun is a function handle to a function 
% % that accepts a vector theta and returns a scalar.
% %
% % Set the third argument to true to regularize the Hessian by adding a
% % multiple of identity. We add sqrt(eps)*eye(size(H))
% %
%     % Set step size.
%     step = eps^(1/4);
%     
%     % Initialize output.
%     p = length(theta);    
%     H = zeros(p,p);
%     for i = 1:p
%         % Use central differences.
%         theta1 = theta;
%         theta1(i) = theta1(i) - step;
%         
%         theta2 = theta;
%         theta2(i) = theta2(i) + step;
%         
%         H(:,i) = ( getGradient(fun, theta2, step) -  getGradient(fun, theta1, step) )/2/step;       
%     end
%     
%     if ( nargin == 3 && wantRegularized )        
%         lambda = eig(H);
%         if ~all(lambda > 0)
%             % There is a -ve or zero eigenvalue.
%             delta = abs(min(lambda)) + sqrt(eps);
%             % Add a multiple of identity to H.
%             H = H + delta*eye(size(H));
%         end
%     end
%     
%     %H = replaceInf(H,realmax);
%     
% end % end of getHessian.

%=== getHessian
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

%=== getDiagonalHessian
function Hdiag = getDiagonalHessian(fun,theta)
% Get a diagonal Hessian of the function fun at theta. We need the diagonal
% Hessian to initialize the B matrix in quasi-Newton optimization. The
% outoput Hdiag is such that:
%
% Hdiag(i,i) = 
% (fun(theta + 2*h*ei) - 2*fun(theta) + fun(theta - 2*h*ei))/(4*h^2)
%
% where ei is a vector of all zeros except a 1 at position i.

    % (1) Set step size.
    step = eps^(1/4);

    % (2) Initialize output.
    p = length(theta);
    Hdiag = zeros(p,p);
    
    % (3) Evaluate fun(theta) once.
    funtheta = fun(theta);
    
    % (4) Get the diagonal elements of Hdiag.
    for i = 1:p                
        % Use central differences.
        theta2 = theta;
        theta2(i) = theta2(i) + 2*step;
        
        theta1 = theta;
        theta1(i) = theta1(i) - 2*step;
        
        Hdiag(i,i) = (fun(theta2) + fun(theta1) - 2*funtheta)/4/step/step;
    end

end % end of getDiagonalHessian.

%=== replaceInf
function B = replaceInf(B,value)
%replaceInf - Replace Inf values in B.
%   B = replaceInf(B,value) takes a matrix or vector B, a scalar value and 
%   replaces +Inf elements in B with abs(value) and -Inf elements in B with
%   -abs(value). The resulting B is returned.

        % (1) Ensure that B is a numeric, matrix and value is a numeric 
        % scalar.
        assert( isnumeric(B) & ismatrix(B) );
        assert( isnumeric(value) & isscalar(value) );

        % (2) Get abs(value).
        absvalue = abs(value);
        
        % (3) Find +Inf or -Inf in B.
        isinfB = isinf(B);
        
        % (4) Find +Inf elements in B and replace them with abs(value).
        B(isinfB & B > 0) = absvalue;
        
        % (5) Find -Inf elements in B and replace them with -abs(value).
        B(isinfB & B < 0) = -absvalue;
end % end of replaceInf.

%=== validateFun
function fun = validateFun(fun)
%validateFun - Validate the fun input to fminqn.
%   fun = validateFun(fun) takes fun which is expected to be a function
%   handle and validates it. If not valid, an error message is thrown.

    % <entry key="BadFun">fun must be a function handle.</entry>
    assertThat(isa(fun,'function_handle'),'stats:classreg:regr:lmeutils:fminqn:BadFun');

end % end of validateFun

%=== validateTheta0
function theta0 = validateTheta0(theta0)
%validateTheta0 - Validate the theta0 input to fminqn.
%   theta0 = validateTheta0(theta0) takes a potential vector theta0 and
%   validates it. If not valid, an error message is thrown. The output
%   vector theta0 will be a column vector.
%
%   What is checked?
%
%   (1) theta0 must be a numeric, real vector.
%
%   (2) theta0 must not contain any NaN or Inf values.

    % <entry key="BadTheta0_NumericRealVector">theta0 must be a numeric, real vector.</entry>
    assertThat(isnumeric(theta0) & isreal(theta0) & isvector(theta0),'stats:classreg:regr:lmeutils:fminqn:BadTheta0_NumericRealVector');
    
    % <entry key="BadTheta0_NoNaNInf">theta0 must not contain any NaN or Inf values.</entry>
    assertThat(~any(isnan(theta0)) & ~any(isinf(theta0)),'stats:classreg:regr:lmeutils:fminqn:BadTheta0_NoNaNInf');
        
    if size(theta0,1) == 1
        theta0 = theta0';
    end

end % end of validateTheta0.

%=== errorIfNotScalar
function errorIfNotScalar(funtheta0)
%errorIfNotScalar - Errors out if funtheta0 is not scalar. Here funtheta0 
% is the function handle fun evaluated at the initial point theta0.
    
    % <entry key="BadFunTheta0">Function handle fun evaluated at theta0 must return a scalar value.</entry>
    assertThat(isscalar(funtheta0),'stats:classreg:regr:lmeutils:fminqn:BadFunTheta0');       

end % end of errorIfNotScalar.

%=== helper method for verifying assertions
function assertThat(condition,msgID,varargin)
%   assertThat(condition,msgID,varargin) takes a variable condition that is
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
