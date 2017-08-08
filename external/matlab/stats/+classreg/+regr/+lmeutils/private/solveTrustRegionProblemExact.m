function [p,mp,isposdef,iter,lambda] = solveTrustRegionProblemExact(g,B,Delta)
%[p,mp,isposdef,iter,lambda] = solveTrustRegionProblemExact(g,B,Delta) 
% computes the exact solution to the "trust region" problem:
% 
%   min w.r.t p : m(p) = g'*p + 0.5* p'*B*p
%
%    subject to : p'*p <= Delta^2 or norm(p,2) <= Delta.
%
%   POSITIONAL PARAMETERS:
%
%           g - A n by 1 real vector.
%
%           B - A n by n real and *symmetric* matrix.
%
%       Delta - The trust region radius, a positive real number. 
%
%   OUTPUTS:
%
%           p - A n by 1 vector which is an estimated solution to the trust 
%               region problem described above.
%
%          mp - Value of the trust region objective m(p) evaluated at p.
%
%    isposdef - Equal to true if B is positive definite and false otherwise.
%
%        iter - Number of iterations of the Newton root finding process 
%               performed (if applicable).
%
%      lambda - The estimated Lagrange multiplier value that satisfies the 
%               necessary and sufficient conditions of optimality.
%
%   
%   NOTES:
%
%    (1) This function assumes that B is *symmetric*. However, we do not check
%        this explicitly in the code. Hence, the responsibility of ensuring 
%        that B is symmetric is on the caller. If B is not symmetric then the
%        results are non-sensical. If B is close to symmetric you can pass
%        in (B+B')/2 instead of B as input to this function.
%
%    (2) This function computes the eigen decomposition of B. Hence, using
%        this function may not be practical if B is a "large" matrix.
%   

%% Problem statement
% We compute the exact solution of the following Trust region subproblem:
%
% $\mbox{min}_{p} \,\,\,\, m(p) = g^T p + \frac{1}{2}  \, p^T B p$
%
% $\mbox{subject to} \,\,\,\, p^T p \le \Delta^2$
%
% Here $B$ is any real symmetric matrix.
%
% We outline the necessary and sufficient conditions of optimality described in
% Nocedal and Wright. 
%
% A vector $p$ is a solution to this subproblem if and only if there exists
% a scalar $\lambda \ge 0$ such that the following conditions are
% satisfied:
%
% (1) $\lambda \ge 0$
%
% (2) $p^T p \le \Delta^2$
%
% (3) $(B + \lambda I) p = -g$
%
% (4) $(B + \lambda I)$ is positive semidefinite
%
% (5) $\lambda (\Delta^2 - p^T p) = 0$ 
%
%
% Consider the following cases:
%
%% Case 1 
% $B$ is symmetric positive semidefinite. Let $g$ be a $n \times 1$ vector and $B$ be a 
% $n \times n$ matrix. Suppose $\lambda_1 \le \lambda_2 \ldots \le
% \lambda_n$ are the real eigenvalues of real symmetric $B$ arranged in an increasing
% order. For Case I, we have $\lambda_1 \ge 0$.
%
%% _Subcase A:_
%
% $\lambda_1 = 0$. Then $B$ is singular.
%
% Let $\hat{p} = - pinv(B) \, g$. 
%
% If $B \hat{p} = -g$ holds with $\hat{p}^T \hat{p} \le \Delta^2$ then $\lambda = 0$ 
% is a possible value of $\lambda$. Otherwise, we conclude that $\lambda >
% 0$. We choose the best bet $\hat{p}$ using the method described in Case 2, Subcase A
% which avoids choosing a threshold in computing pinv.
%
% If $\lambda > 0$ then $(B + \lambda I)$ is positive definite. The step
% length is given by:
%
% $p(\lambda) = -(B + \lambda I)^{-1} g$ with $\lambda$ given by:
%
% $f(\lambda) = p(\lambda)^T p(\lambda) = \Delta^2$.
% 
% Our plan is to solve this equation for $\lambda$ by Newton's method. It
% is easier to solve a different version of the above equation:
%
% $\frac{1}{ \sqrt{ f(\lambda) } } = \frac{1}{ \Delta }$.
% 
% Since we want to keep $\lambda > 0$ during Newton iterations we
% parameterize $\lambda = e^{\theta}$. Suppose $B = Q \Lambda Q^T$ is the
% eigen decomposition of $B$ where $Q = [q_1,q_2,\ldots,q_n]$ and $\Lambda
% = diag(\lambda_1,\lambda_2,\ldots,\lambda_n)$. $q_i$ is an eigenvector
% corresponding to $\lambda_i$. It is easy to see that:
%
% $p(\lambda) = -\sum_{i=1}^n q_i \frac{(q_i^T g)}{(\lambda_i + \lambda)}$
%
% and 
%
% $f(\lambda) = p(\lambda)^T p(\lambda) = \sum_{i=1}^n \frac{ (q_i^T g)^2
% }{ (\lambda_i + \lambda)^2 }$
%
% and
%
% $f(\theta) = \sum_{i=1}^n \frac{ (q_i^T g)^2 }{ (\lambda_i + e^{\theta})^2 }$
%
% Let $\phi(\theta) = \frac{1}{ \sqrt{ f(\theta) } } - \frac{1}{ \Delta }$.
% 
% We would like to find zeros of $\phi(\theta)$ using Newton's method.
%
%
% $\phi'(\theta) = -\frac{1}{2} \frac{ f'(\theta) }{ f(\theta) f(\theta)^{1/2} }$
%
% The Newton update is given by:
%
% $\theta_{new} = \theta - \frac{\phi(\theta)}{\phi'(\theta)}$
%
% Substituting and simplifying:
%
% $\theta_{new} = \theta + \frac{\Delta - \sqrt{f(\theta)}}{\Delta} \left\{ \frac{2 f(\theta)}{f'(\theta)} \right\}$
%
% where
%
% $f'(\theta) = (-2 e^{\theta}) \sum_{i=1}^n \frac{(q_i^T g)^2}{(\lambda_i
% + e^{\theta})^3} = (-2 e^{\theta}) \,\, q(\lambda)^T q(\lambda) = (-2 e^{\theta}) \,\, p(\lambda)^T (B + \lambda I)^{-1} p(\lambda)$
%
% with $\lambda = e^{\theta}$.
% 
% If $(B + \lambda I) = L L^T$ then we can write:
%
% $L L^T p(\lambda) = -g$ and
%
% $L q(\lambda) = p(\lambda)$ with
%
% $f'(\theta) = (-2 e^{\theta}) \,\, q(\lambda)^T q(\lambda)$.
%
% Substituting $f(\theta)$ and $f'(\theta)$ in the Newton update equation we get: 
%
% $\theta_{new} = \theta + \frac{norm(p(\lambda),2) - \Delta}{ (e^{\theta}) \Delta} \left\{ \frac{ p(\lambda)^T p(\lambda)}{ q(\lambda)^T q(\lambda)} \right\}$
%
%% _Subcase B:_
%
% If $\lambda_1 > 0$. Then $B$ is positive definite.
%
% Let $\hat{p} = B^{-1} g$. If $\hat{p}^T \hat{p} \le \Delta^2$ then
% $\lambda = 0$ will work. Otherwise, $\lambda > 0$ and we proceed as for
% the case $\lambda > 0$ in Subcase A of Case 1.
%
%
%% Case 2
% $B$ is symmetric but not positive definite. 
%
% The smallest eigenvalue of $B$ satisfies $\lambda_1 < 0$. Since $B +
% \lambda I$ must be positive semidefinite, we must have $\lambda \ge
% -\lambda_1 > 0$.
%
%% _Subcase A:_
%
% If $\lambda = -\lambda_1$ then $g$ must satisfy $q_j^T g = 0$ for all $j$ with
% $\lambda_j = \lambda_1$ in order for $(B -\lambda_1 I) p = -g$ to have a solution. 
% Let $\hat{p} = -pinv( B -\lambda_1 I ) g$. Assuming $(B -\lambda_1 I) \hat{p} =
% -g$ holds, one general solution is given by:
%
% $p = \hat{p} + \tau \, q_1$
%
% where
%
% $\hat{p} = -\sum_{i=1, \lambda_i \neq \lambda_1}^n q_i \frac{(q_i^T g)}{(\lambda_i -\lambda_1)}$
%
% Also
%
% $p^T p = \hat{p}^T \hat{p} + \tau^2$ since $\hat{p}$ and $q_1$ are
% orthogonal.
%
% Since $\lambda > 0$ we must have $p^T p = \Delta^2$. This gives us a
% value for $\tau$ and the solution $p$.
%
% In general, if $\lambda_1$ has multiplicity $k$ (i.e., $\lambda_2 = \ldots = \lambda_k = \lambda_1$)
% then the solution $p$ is given by:
%
% $p = \tau \, q_1 + \sum_{i=k+1}^n q_i \frac{(-q_i^T g)}{(\lambda_i -\lambda_1)}$
%
% and so
%
% $p^T p = \tau^2 + \sum_{i=k+1}^n \frac{(q_i^T g)^2}{(\lambda_i -\lambda_1)^2} = \Delta^2$
%
% If $\lambda = -\lambda_1$ holds then our best bet for $p$ is given by the
% expression above with the multiplicity $k$ chosen to be the smallest
% value that satisfies:
%
% $\sum_{i=k+1}^n \frac{(q_i^T g)^2}{(\lambda_i -\lambda_1)^2} \le \Delta^2$
%
% Why the smallest $k$? From the solution for $p$ it is clear that:
%
% $q_m^T p \, (\lambda_m -\lambda_1) = - q_m^T g$ for $m \ge (k+1)$.
%
% If $(B -\lambda_1 I) p = -g$ holds then we will have:
%
% $q_m^T p \, (\lambda_m -\lambda_1) = - q_m^T g$ for $m = 1,\ldots,n$.
%
% Hence by choosing $k$ to be as small as possible, we try to satisfy the
% set of $n$ equations above as closely as possible.
%
% Once we determine the best bet for $p$, we check if $(B -\lambda_1 I) p = -g$ 
% holds or not.
%
% Having determined smallest $k$ such that $\sum_{i=k+1}^n \frac{(q_i^T
% g)^2}{(\lambda_i -\lambda_1)^2} \le \Delta^2$, we compute:
%
% $\tau_1 = \sqrt{ \Delta^2 - \sum_{i=k+1}^n \frac{(q_i^T g)^2}{(\lambda_i
% -\lambda_1)^2} }$
%
% $\tau_2 = - \sqrt{ \Delta^2 - \sum_{i=k+1}^n \frac{(q_i^T g)^2}{(\lambda_i
% -\lambda_1)^2} }$
%
% There are many possible solutions for $p$ including:
%
% $p = \tau_1 \, q_1 + \sum_{i=k+1}^n q_i \frac{(-q_i^T g)}{(\lambda_i -\lambda_1)}$
%
% $p = \tau_2 \, q_1 + \sum_{i=k+1}^n q_i \frac{(-q_i^T g)}{(\lambda_i -\lambda_1)}$
%
% $p = \tau_1 \, q_2 + \sum_{i=k+1}^n q_i \frac{(-q_i^T g)}{(\lambda_i -\lambda_1)}$
%
% $p = \tau_2 \, q_2 + \sum_{i=k+1}^n q_i \frac{(-q_i^T g)}{(\lambda_i -\lambda_1)}$
%
% $\vdots$
%
% $p = \tau_1 \, q_k + \sum_{i=k+1}^n q_i \frac{(-q_i^T g)}{(\lambda_i -\lambda_1)}$
%
% $p = \tau_2 \, q_k + \sum_{i=k+1}^n q_i \frac{(-q_i^T g)}{(\lambda_i -\lambda_1)}$
%
% Theoretically, these should give the same value for the objective function $m(p)$.
% However, in floating point arithmetic, this is not so. Hence, we compute
% $m(p)$ for each potential $p$ and choose the $p$ that gives the smallest
% $m(p)$.
%
% If $\lambda = -\lambda_1 = 0$ then $p^T p < \Delta^2$ is acceptable. In
% this case, we must also evaluate the candidate:
%
% $p = \sum_{i=k+1}^n q_i \frac{(-q_i^T g)}{(\lambda_i -\lambda_1)}$
%
% In the above discussion, $k = 0$ is acceptable corresponding to
% multiplicity $0$.
%
%% _Subcase B:_
%
% If $\lambda > -\lambda_1$ then $B + \lambda I$ is positive definite. We
% can solve:
%
% $p(\lambda)^T p(\lambda) = \Delta^2$ while ensuring $\lambda >
% -\lambda_1$. We can do this by parameterizing:
%
% $\lambda = -\lambda_1 + e^{\theta}$.
%
% The Newton update equation is the same as before:
%
% $\theta_{new} = \theta + \frac{norm(p(\lambda),2) - \Delta}{ (e^{\theta}) \Delta} \left\{ \frac{ p(\lambda)^T p(\lambda)}{ q(\lambda)^T q(\lambda)} \right\}$
%
%
%% Choosing initial theta0
%
% How to choose initial $\theta_0$ for starting Newton Iterations? We
% present a simple heuristic that appears to work well. First, we
% approximate $B$ using a diagonal matrix $\delta I$. We choose $\delta$
% so as to minimize the Frobenius norm of the difference between $B$ and $\delta I$:
%
% $\delta = \mbox{arg min}_{\delta}  \,\, \mbox{trace}\{ (B-\delta I)^T (B-\delta I) \}$
%
% The solution is given by:
%
% $\delta = \mbox{trace}(B)/n$
%
% With the approximation, $B = \delta I$, we would like to choose $p$ and $\lambda_0$ to
% satisfy $(B + \lambda_0 I)p = -g$ such that $p^T p = \Delta^2$.
%
% First, note that $p$ is simply given by:
%
% $p = -g/(\delta + \lambda_0)$
%
% and with this choice $p^T p = \Delta^2$ implies:
%
% $g^T g/(\delta + \lambda_0)^2 = \Delta^2$
%
% Solving for $\lambda_0$ we get:
%
% $\lambda_0 = -\delta + norm(g)/\Delta$ or
%
% $\lambda_0 = -\delta - norm(g)/\Delta$
%
% However, the requirement that $(B + \lambda_0 I) = (\delta + \lambda_0)
% I$ be positive semidefinite means that $\lambda_0 \ge -\delta$. Hence the
% only acceptable solution for $\lambda_0$ is:
%
% $\lambda_0 = -\delta + norm(g)/\Delta$
%
% There are 2 other constraints that $\lambda_0$ must satisfy:
%
% $\lambda_0 \ge 0$ and
% 
% $\lambda_0 \ge -\lambda_1$
%
% Hence, I propose to choose $\lambda_0$ as follows:
%
% $\lambda_0 = \mbox{max}\{ -\delta + norm(g)/\Delta, \,\,\, 0, \,\,\, -\lambda_1 \}$
%
% If the true value of $\lambda$ is $0$, we wouldn't want the Newton
% iteration to start at exactly $0$ but rather at some small value.
% Similarly, if the true value of $\lambda$ is $-\lambda_1$ we would like
% the Newton iteration to start at a value close to $-\lambda_1$ but not
% exactly at $-\lambda_1$. Hence, we define our choice of $\lambda_0$ as
% follows:
%
% $\lambda_0 = \mbox{max}\{ -\delta + norm(g)/\Delta, \,\,\, \gamma \, 100 \, \sqrt{absTol}, \,\,\, -\gamma \, \lambda_1 \}$
%
% where $absTol$ is the absolute tolerance for testing against $0$ and
% $\gamma$ is a factor close to 1.0. I choose $\gamma = 1.2$ but other
% choices work equally well.
%
% Having chosen $\lambda_0$, the initial $\theta_0$ using the relation $\lambda = -\lambda_1 + e^{\theta}$
% is given by:
%
% $\theta_0 = \log( \lambda_0 + \lambda_1 )$
%

%   Copyright 2012-2013 The MathWorks, Inc.    

%% solveTrustRegionProblemExact   

    % Compute the eigen decomposition of n by n matrix B, sort the eigenvalues 
    % in ascending order and rearrange the corresponding eigenvectors. After 
    % this process, we will get a n by 1 column vector Lambda and a n by n 
    % matrix Q such that:
    %   Lambda(i) = ith largest eigenvalue and 
    %      Q(:,i) = eigenvector corresponding to Lambda(i).
         [Q, Lambda] = eig(B);
              Lambda = diag(Lambda);
              
         % Sort only if you need to.     
         if ~issorted(Lambda)
            [Lambda, sidx] = sort(Lambda,'ascend');
                         Q = Q(:,sidx);
         end
                   
    % Precompute certain quantities for later use.
    %
        % Absolute values of eigenvalues.              
        absLambda = abs(Lambda);
       
        % Largest eigenvalue in absolute value.
        maxabsLambda = max(absLambda);
   
        % Store the product Q'*g for later use.
        Qtg = Q'*g;
    
    % Convergence tolerances.    
    %    
        % Relative error tolerance when testing a quantity against zero. 
        % For double data, this is around 1.8190e-12.
        relTol = eps(class(B))^(3/4);

        % Absolute error tolerance when testing a quantity against zero.
        % For double data, this is around 1.8190e-12.
        absTol = eps(class(B))^(3/4);
        
        % In some cases, we will be required to find a step length p such that
        % norm(p) = Delta. We will assume that the required p has been found if
        % abs( norm(p)/Delta - 1 ) <= trustTol. For double data, this is
        % around 1.4901e-08.
        trustTol = eps(class(B))^(1/2);
        
    % Process the various possibilities based on Lambda(1).
    %
        if ( (Lambda(1) > 0 && abs(Lambda(1)) <= relTol * maxabsLambda) || Lambda(1) == 0 ) 
            % Case 1. Subcase A. Lambda(1) = 0.

                % Mark B as not positive definite.
                isposdef = false;
                                                                    
                % (1) Assume lambda = 0 does not hold, then lambda > 0.
                % (B + lambda*I) is positive definite for lambda > 0.
                % Parameterize: lambda = 0 + exp(theta).
                       offset = 0;
                    LambdaFun = makeLambdaFun( offset );
                    
                    % Starting value of theta.
                    theta0 = chooseInitialTheta(g,B,Delta,Lambda(1),absTol);

                    % Find a value of theta such that p = inv(B + lambda*I) * (-g)
                    % satisfies ||p||_2 = Delta where lambda = LambdaFun( theta ).                         
                    [theta,p,iter,cause] = solveStepLengthEquation(Qtg, Q, Lambda, Delta, LambdaFun, theta0, trustTol);
                   
                    % Compute lambda corresponding to this theta.
                    lambda = LambdaFun(theta);                                      
 
                % (2) Try lambda = 0. If root finding did not converge (cause ~= 0)
                % then p may contain some NaN or Inf values. Assume Subcase A 
                % holds with lambda = 0 and re-compute p if:
                %   (a) lambda from root-finding is "small". 
                %                      or
                %   (b) Delta is Inf.
                %                      or
                %   (c) norm(p) == 0. This happens when all elements of Qtg
                %   are 0. We could as well check for all(Qtg == 0) instead
                %   of norm(p) == 0 for clarity.
                
                    if ( cause ~= 0 )
                        % Root finding did not converge.
                        
                        % Is lambda from root finding small enough?
                        islambdaSmall = abs(lambda) <= sqrt(absTol);
                        
                        % Compute p again assuming lambda = 0 if required.     
                            recomputep = islambdaSmall || isinf(Delta) || (norm(p) == 0);                        
                            if ( recomputep )
                                % norm(p) = Delta not required when lambda = 0.
                                     p = bestBetForStep(g, B, Qtg, Q, Lambda, 0, Delta);
                                lambda = 0;
                            else
                                % Warn that the problem was not solved.
                                    warning(message('stats:classreg:regr:lmeutils:fminqn:Message_UnableToSolveTrustRegionProblem'));
                                    %warning('solveTrustRegionProblemExact:ProblemNotSolved','Unable to solve problem.');
                                     p = NaN(length(g),1);                                                                        
                            end
                    end
                                    
                % (3) Compute trust region objective at p.
                mp = trustRegionObjective(g,B,p);                    
                   
        elseif ( Lambda(1) > 0 && abs(Lambda(1)) > relTol * maxabsLambda )
            % Case 1. Subcase B. Lambda(1) > 0.

                % Mark B as positive definite.
                isposdef = true;
            
                % (1) Try lambda = 0. 
                    [phat, kbst] = bestBetForStep(g, B, Qtg, Q, Lambda, 0, Delta);
                    if ( kbst == 0 )
                        % If kbst = 0, then lambda = 0 works
                        % with B*p = -g and norm(p) <= Delta.
                               p = phat;
                              mp = trustRegionObjective(g,B,p);
                            iter = 0;
                          lambda = 0;
                          return;
                    end
                    
                % (2) Assume lambda = 0 does not hold, then lambda > 0.
                % Since Lambda(1) > 0, B is positive definite and so 
                % (B + lambda*I) should also be positive definite for 
                % lambda >= 0. Parameterize lambda = 0 + exp(theta).
                       offset = 0;
                    LambdaFun = makeLambdaFun( offset );

                    % Starting value of theta.
                    theta0 = chooseInitialTheta(g,B,Delta,Lambda(1),absTol);

                    % Find a value of theta such that p = inv(B + lambda*I) * (-g)
                    % satisfies ||p||_2 = Delta where lambda = LambdaFun( theta ).                             
                    [theta,p,iter,cause] = solveStepLengthEquation(Qtg, Q, Lambda, Delta, LambdaFun, theta0, trustTol);

                    % Compute lambda corresponding to this theta.
                    lambda = LambdaFun(theta);  
                
                % (3) Try lambda = 0. If root finding did not converge (cause ~= 0)
                % then p may contain some NaN or Inf values. Assume Subcase B
                % holds with lambda = 0 and re-compute p if:
                %   (a) lambda from root-finding is "small". 
                %                     or
                %   (b) Delta is Inf.
                %                     or
                %   (c) norm(p) == 0. This happens when all elements of Qtg
                %   are 0. We could as well check for all(Qtg == 0) instead
                %   of norm(p) == 0 for clarity.
                
                    if ( cause ~= 0 )
                        % Root finding did not converge.
                        
                        % Is lambda from root finding small enough?
                        islambdaSmall = abs(lambda) <= sqrt(absTol);
                        
                        % Compute p again assuming lambda = 0 if required.     
                            recomputep = islambdaSmall || isinf(Delta) || (norm(p) == 0);                        
                            if ( recomputep )
                                % norm(p) = Delta not required when lambda = 0.
                                % We are basically doing this: p = -1 * Q * ( Qtg ./ Lambda );
                                     p =  phat;
                                lambda = 0;
                            else
                                % Warn that the problem was not solved.
                                    warning(message('stats:classreg:regr:lmeutils:fminqn:Message_UnableToSolveTrustRegionProblem'));
                                    %warning('solveTrustRegionProblemExact:ProblemNotSolved','Unable to solve problem.');
                                     p = NaN(length(g),1);                                                                       
                            end
                    end                   
                    
                % (4) Compute trust region objective at p.
                mp = trustRegionObjective(g,B,p);

        elseif ( Lambda(1) < 0 )
            % Case 2. Lambda(1) < 0.

                % Mark B as not positive definite.
                isposdef = false;
                            
                % (1) Subcase B. Assume lambda = -Lambda(1) does not hold. 
                % Then lambda > -Lambda(1).
                %
                    % (B + lambda*I) is positive definite for lambda > -Lambda(1).
                    % Parameterize lambda = -Lambda(1) + exp(theta).
                       offset = -Lambda(1);
                    LambdaFun = makeLambdaFun( offset );

                    % Starting value of theta.
                    theta0 = chooseInitialTheta(g,B,Delta,Lambda(1),absTol);

                    % Find a value of theta such that p = inv(B + lambda*I) * (-g)
                    % satisfies ||p||_2 = Delta where lambda = LambdaFun( theta ).       
                    [theta,p,iter,cause] = solveStepLengthEquation(Qtg, Q, Lambda, Delta, LambdaFun, theta0, trustTol);

                    % Compute lambda corresponding to this theta.
                    lambda = LambdaFun(theta);    

                % (2) Subcase A. Try lambda = -Lambda(1). If root finding did 
                % not converge (cause ~= 0) then p may contain some NaN or Inf values. 
                % Assume lambda = -Lambda(1) and re-compute p if:
                %   (a) lambda from root-finding is close to -Lambda(1). 
                %                   or
                %   (b) norm(p) == 0. This happens when all elements of Qtg
                %   are 0. We could as well check for all(Qtg == 0) instead
                %   of norm(p) == 0 for clarity.

                    if ( cause ~= 0 )
                        % Root finding did not converge.
                        
                         % Is lambda from root finding close to -Lambda(1)?
                         islambdaPlusLambda1Small = abs(lambda + Lambda(1)) <= sqrt(absTol)*max(1,abs(Lambda(1)));
                        
                         % Compute p again assuming lambda = -Lambda(1) if required.     
                         recomputep = islambdaPlusLambda1Small || (norm(p) == 0);                                  
                         if ( recomputep )
                            % norm(p) = Delta is required when lambda = -Lambda(1) > 0.
                                 p = bestBetForStep(g, B, Qtg, Q, Lambda, -Lambda(1), Delta);
                            lambda = -Lambda(1);
                         else
                            % Warn that the problem was not solved.
                                warning(message('stats:classreg:regr:lmeutils:fminqn:Message_UnableToSolveTrustRegionProblem'));
                                %warning('solveTrustRegionProblemExact:ProblemNotSolved','Unable to solve problem.');
                                 p = NaN(length(g),1);                                 
                         end                        
                    end                                       
                    
                % (3) Compute trust region objective at p.
                mp = trustRegionObjective(g,B,p);
                                                
        end % end of if.                  
        
        % Final sanity check. mp must be negative.
        if mp > 0 || isnan(mp)
            p = zeros(length(g),1);
            mp = 0;
        end        
        
end % end of solveTrustRegionProblemExact.

%% bestBetForStep.
%=== bestBetForStep
function [p,kbst,mbest] = bestBetForStep(g, B, Qtg, Q, Lambda, lambda, Delta)
% Find a best bet solution p to the problem:
%
% (1) (B + lambda*I)*p = -g such that
% (2)             p'*p = Delta^2
%
% => The eigenvalues of B are in the vector Lambda arranged in ascending order. 
% => The columns of Q contain the eigenvectors corresponding to eigenvalues in Lambda. 
% => We are given the vector Q'*g. 
% => We are also given the values of lambda and Delta.
%
% Our assumption is that the input value lambda is such that:  
% Lambda(1) + lambda, Lambda(2) + lambda, ..., Lambda(k) + lambda are effectively zero.
% In other words, Lambda(1) has effective multiplicity k.
%
% We determine the multiplicity k and a best bet solution p that satisfies 
% eqn. (2) exactly and also 
% 
%    Q(:,L)^'*p = -Q(:,L)'*g/(Lambda(L) + lambda) holds for L >= (k+1) by construction.
%
% The only thing we are not sure of is whether Q(:,M)'*g = 0 for M <= k holds or not.
% But this can be easily checked by checking if (B + lambda*I)*p = -g holds. If so, 
% then p satisfies both (1) and (2).
% 

    % (1) Let x = (Qtg).^2 ./ (Lambda + lambda).^2.
    % Find the smallest value of k such that sum( x(k+1:end) ) <= Delta^2. 
    % For this value of k determine tau2 = Delta^2 - sum( x(k+1:end) ).
        n = length(Qtg);
        %x = ( Qtg.^2 ) ./ ( (Lambda + lambda).^2 ); 
        x = ( Qtg./(Lambda + lambda) ).^2; 
        for k = 0:n
            kbst = k;
            tau2 = Delta^2 - sum( x(k+1:end) );            
            if ( tau2 >= 0 )
                % We found k (in kbst) and the corresponding tau2.
                break;
            end
        end
    
     % (2) Compute p = c1 + c2 where c2 is the fixed component of p and c1 
     % is a variable component of p.        
     %   
     
        % The fixed component of p.
            idx = (kbst+1):n;
            if ~isempty(idx)
                c2 = -1 * Q(:,idx) * ( Qtg(idx) ./ (Lambda(idx) + lambda) );
            else
                c2 = zeros(n,1);
            end
            
        % If kbst == 0, c2 is the required solution.
            if ( kbst == 0 )
                    p = c2;
                mbest = trustRegionObjective(g,B,p);
                return;
            end
            
        % The variable component c1 of p looks like:
        %
        %  c1 = sqrt(tau2)*Q(:,m) or -sqrt(tau2)*Q(:,m) for m = 1..kbst.
        %
        % There are a total of 2*kbst possibilities for c1. We choose the
        % variable component c1 such that trustRegionObjective(g,B,c1+c2)
        % is minimized. 
        
            
            % Best step and objective found so far. If Delta or tau2 is
            % Inf, we do not want to force selection of c2. Hence make
            % pbest and mbest NaN and Inf respectively in this case.
            if ( lambda == 0 && ~isinf(Delta) && ~isinf(tau2) )
                % We can have p'*p < Delta^2. Hence c2 is also a
                % possibility for p.
                pbest = c2;
                mbest = trustRegionObjective(g,B,pbest);   
            else
                pbest = NaN(n,1);
                mbest = Inf;
            end
                
            % Compare various possibilities for p and choose the best one.    
                for m = 1:kbst       
                    % Choose best possibility for current m.
                    %
                        % Two possible values of c1: c1pos and c1neg.
                            c1pos =  sqrt(tau2)*Q(:,m);
                            c1neg = -sqrt(tau2)*Q(:,m);

                        % Two possible values of p: ppos and pneg.
                            ppos = c1pos + c2;
                            pneg = c1neg + c2;

                        % Possible values of trust region objective: mpos and mneg.
                            mpos = trustRegionObjective(g,B,ppos);
                            mneg = trustRegionObjective(g,B,pneg);

                        % Set pcurr and mcurr by comparing (ppos,mpos) against (pneg,mneg).    
                            if ( mpos <= mneg )
                                pcurr = ppos;
                                mcurr = mpos;
                            else
                                pcurr = pneg;
                                mcurr = mneg;
                            end

                    % Compare (pcurr,mcurr) with (pbest,mbest).
                        if ( mcurr <= mbest )
                            mbest = mcurr;
                            pbest = pcurr;
                        end
                end % end of for.
        
            % Set p equal to pbest.
            p = pbest;
     
end % end of bestBetForStep.

%% solveStepLengthEquation.
%=== solveStepLengthEquation
function [theta,p,iter,cause] = solveStepLengthEquation(Qtg, Q, Lambda, Delta, LambdaFun, theta0, trustTol)
% Find a value of theta such that p = inv(B + lambda*I) * (-g)
% satisfies ||p||_2 = Delta where lambda = LambdaFun( theta ).
%
% INPUTS:
%        Qtg = n by 1 vector containing the product Q'*g.
%          Q = n by n matrix where each column contains an eigenvector of B.
%     Lambda = n by 1 vector containing sorted eigenvalues of B. Lambda(i) corresponds to Q(:,i). 
%      Delta = real positive trust region radius.
%  LambdaFun = function handle to a function that outputs lambda given an
%              input value theta.
%     theta0 = initial value of theta.
%   trustTol = relative convergence threshold for stopping iterations.
%
% OUTPUTS:
%      theta = optimal value of theta.
%          p = optimal value of p.
%       iter = number of iterations taken for convergence.
%      cause = reason for termination.
%
% cause is indicated by integer codes as described below:
%
%   cause = 0 -> abs( sqrt(ptp)/Delta - 1 ) <= trustTol was satisfied.
%   cause = 1 -> abs(lambdaNew - lambda)    <= trustTol was satisfied.
%   cause = 2 -> iter > maxit, iteration limit reached.
%   cause = 3 -> isnan(NewtonStep) || isinf(NewtonStep), NaN or Inf NewtonStep.

    % Convergence flag.
    found = false;   
    
    % Initialize theta, lambda, p'*p and zeroViolation. We will always
    % keep these 4 quantities synchronized.
    %
        % Current value of theta.
        theta = theta0;
    
        % Current value of lambda.
        lambda = LambdaFun(theta);                
        
        % Compute p'*p.
        %ptp = sum( (Qtg.^2) ./ ((Lambda + lambda).^2) );
        ptp = sum( ( Qtg./(Lambda + lambda) ).^2 );
    
        % Current difference between 1/sqrt(ptp) and 1/Delta. We will call
        % this the zero violation.
        zeroViolation = abs( 1/sqrt(ptp) - 1/Delta );
    
    % Iteration counter and maximum number of iterations allowed.
        iter = 0;
       maxit = 100;
    
    % Address trivial case.
        if ( Delta == 0 )
                p = zeros(length(Qtg),1);
            cause = 0;
            return;
        end       
       
    % Start Newton iterations.
    while( found == false )
        
            % Compute qtq.
            %qtq = sum( (Qtg.^2) ./ ((Lambda + lambda).^3) );
            qtq = sum( (Qtg./(Lambda + lambda)).^2 ./ (Lambda + lambda) );
            
            % Compute Newton step.
            NewtonStep = ((sqrt(ptp) - Delta)/Delta)*(exp(-theta)*ptp)/qtq;                                    
            
            if ( isnan(NewtonStep) || isinf(NewtonStep) )
                % Terminate with cause = 3 if NewtonStep is not valid.
                found = true;
                cause = 3;
            else
                % NewtonStep is OK.
                
                % Try various stepLength values starting from 1.0. If we take
                % the step then zeroViolation after taking the step must decrease, 
                % otherwise we decrease the stepLength.
                    stepLength = 1.0;
                    acceptStep = false;            
                    while( acceptStep == false )

                            % Take the step and compute new values of theta,
                            % lambda, ptp and zeroViolation.
                                    thetaNew = theta + stepLength * NewtonStep;
                                   lambdaNew = LambdaFun(thetaNew);
                                      %ptpNew = sum( (Qtg.^2)./ ((Lambda + lambdaNew).^2) );
                                      ptpNew = sum( ( Qtg./(Lambda + lambdaNew) ).^2 );
                            zeroViolationNew = abs( 1/sqrt(ptpNew) - 1/Delta );   

                            % Take the step if zeroViolation decreases and ptpNew is sensible.
                                if ( zeroViolationNew <= zeroViolation && ~isinf(ptpNew) )                            
                                    % Mark current step as accepted and update
                                    % theta, lambda, ptp and zeroViolation.
                                       acceptStep = true;                  
                                            theta = thetaNew;
                                           lambda = lambdaNew;
                                              ptp = ptpNew;
                                    zeroViolation = zeroViolationNew;                                                                        
                                else                            
                                    % Reduce stepLength by half and try again.
                                    stepLength = stepLength/2;

                                    % Terminate if lambda is not changing significantly.
                                    if ( abs(lambdaNew - lambda) <= trustTol )
                                       acceptStep = true;
                                            found = true;
                                            cause = 1;
                                    end
                                end

                    end % end of while.
            
            end % end of if.
                
            % Check convergence. Stop iterations if norm(p) = Delta holds
            % or if we exceed the allowed number of iterations.               
                if ( iter > maxit )
                    found = true;
                    cause = 2;
                end                      
            
                if ( abs( sqrt(ptp)/Delta - 1 ) <= trustTol )
                   found = true;
                   cause = 0;
                end           
                
            % Update iteration counter.
            iter = iter + 1;
            
    end % end of while.
    
    % Compute final step size using final lambda.
    p = -1 * Q * ( Qtg./(Lambda + lambda) );        
        
end % end of solveStepLengthEquation.

%% chooseInitialTheta
%=== chooseInitialTheta
function theta0 = chooseInitialTheta(g,B,Delta,Lambda1,absTol)
% Choose initial theta0 for starting Newton iterations.
% INPUTS:
%       g - A n by 1 real vector.
%       B - A n by n real and *symmetric* matrix.
%   Delta - The trust region radius, a positive real number. 
% Lambda1 - The smallest eigenvalue of B.
%  absTol - Absolute tolerance for testing against 0.
%
% OUTPUTS: 
% The initial value theta0.

    % Assuming B is diagonal, the best approximation to B 
    % is given by delta*I where delta = trace(B)/n where n
    % is the size of B. 
            n = length(g);
        delta = trace(B)/n;
    
    % If delta*I is used as an approximation to B then the requirement that
    % (delta*I + lambda0*I) be positive semidefinite implies that lambda0 >= -delta.
    % Using the constraint p'*p = Delta^2, we get:
       lambda0 = -delta + (norm(g)/Delta);

    % Since lambda0 must also be >= 0, >= -Lambda1, we use:
         gamma = 1 + 0.2;
           tol = sqrt(absTol)*100;
       lambda0 = max( [lambda0, gamma*tol, -gamma*Lambda1] );

    % Since lambda = -Lambda1 + exp(theta), so theta = log(lambda + Lambda1).
    theta0 = log(lambda0 + Lambda1);
        
end % end of chooseInitialTheta.

%% makeLambdaFun.
%=== makeLambdaFun
function LambdaFun = makeLambdaFun( offset )
% Return a function handle to a function that accepts a scalar theta and
% returns offset + exp(theta).

    LambdaFun = @f;

        function val = f( theta )
        
            val = offset + exp(theta);
            
        end
    
end % end of makeLambdaFun.

%% trustRegionObjective.
%=== trustRegionObjective
function mp = trustRegionObjective(g,B,p)
% Given g, B and p, compute m(p) = g'*p + 0.5* p'*B*p

    mp = g'*p + 0.5* (p'*B*p);

end % end of trustRegionObjective.
