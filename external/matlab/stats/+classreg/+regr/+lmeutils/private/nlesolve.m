function [xstar,rxstar,Jxstar,cause] = nlesolve(rfun,x0,varargin)
% NLESOLVE - Solves a system of non linear equations.
%   [xstar,rxstar,Jxstar] = nlesolve(rfun,x0) takes a function handle rfun,
%   q-by-1 initial point x0 and solves the system of non linear equations
%   specified by rfun. The function handle rfun should be callable like
%   [rx,Jx,pnx] = rfun(x) where rx is the q-by-1 vector of non linear
%   equations evaluated at q-by-1 vector x and Jx is q-by-q Jacobian of rx
%   evaluated at x. pnx is the solution to (Jx'*Jx)*pnx = -Jx'*rx. Output
%   xstar is such that [rxstar,Jxstar] = rfun(xstar) and either rxstar is
%   approximately 0 or Jxstar'*rxstar is approximately 0. The output cause
%   is an integer indicating reason for termination:
%
%   cause                         Meaning
%     0                       Relative infinity norm of Jxstar'*rxstar is below TolFun.
%     1                       Two norm of step size is below TolX.
%     2                       Iteration limit reached.
%
%   cause = 0 or 1 indicates success and cause = 2 indicates failure.
%
%   [xstar,rxstar,Jxstar] = nlesolve(rfun,x0,varargin) accepts optional
%   name/value pairs as follows:
%
%       'Options'   A structure containing optimization options with the 
%                   following fields:
% 
%            TolFun   - Relative tolerance on the gradient of the objective
%                       function. Default is 1e-6.
%                    
%              TolX   - Absolute tolerance on the step size. Default is 
%                       1e-12.
%  
%           MaxIter   - Maximum number of iterations allowed. Default is
%                       1000.
%  
%           Display   - Level of display:  'off', 'iter', or 'final'.
%                       Default is 'off'. Both 'iter' and 'final' have the
%                       same effect.
%
%       'Method'    A string selecting the solution method. Choices are
%                   'TrustRegion2D' (default), 'LineSearchModifiedNewton',
%                   and 'LineSearchNewton'. 'TrustRegion2D' uses a trust
%                   region method with 2D subspace minimization,
%                   'LineSearchModifiedNewton' uses a line search method
%                   with modified Newton direction and 'LineSearchNewton'
%                   uses a line search method with Newton direction.
%                   'LineSearchNewton' may fail if the Jacobian of the
%                   system of non linear equations becomes singular during
%                   iterations.

%% Problem description
% Suppose $x \in R^q$ and $r: R^q \rightarrow R^q$ be a system of $q$ non
% linear equations. We wish to find $x^*$ such that:
%%
%
% $$r(x^*) = 0$$
%
% We attempt to solve this problem by considering the problem of minimizing
% the following function:
%
% $f(x) = \frac{1}{2} r(x)^T r(x)$
%
% If $J(x)$ is the Jacobian of $r(x)$, the gradient of $f(x)$ is given by:
%
% $\nabla f(x) = J(x)^T r(x)$
%
% If $x^*$ minimizes $f(x)$ then: 
%
% $\nabla f(x^*) = J(x^*)^T r(x^*) = 0$
%
% If $J(x^*)$ is non singular then this implies that $r(x^*) = 0$.
% Therefore solving for $r(x^*) = 0$ by minimizing $f(x)$ is a sensible
% approach as long as $J(x^*)$ is non singular.
%
%% Solution approach
% Our solution approach consists of minimizing $f(x)$ using a trust region 
% technique. We construct a quadratic approximation to $f(x)$ around an
% iterate $x_k$ that we "trust" in a trust region of radius $\Delta_k$. Our
% quadratic approximation is:
%
% $f(x_k + p) = \frac{1}{2} r(x_k + p)^T r(x_k + p)$
%
% $r(x_k + p) \approx r(x_k) + J(x_k) p$
%
% $$m(x_k + p) = \frac{1}{2} \{r(x_k) + J(x_k) p\}^T \{r(x_k) + J(x_k) p\}$
%
% Simplifying, we can write:
%
% $$m(x_k + p) = f(x_k) + r(x_k)^T J(x_k) p + \frac{1}{2} p^T J(x_k)^T J(x_k) p$
% 
% We trust $m(x_k + p)$ when $x_k + p$ is in a radius $\Delta_k$ around
% $x_k$. This implies that $p^T p \le \Delta_k^2$. Since we are looking to
% minimize $f(x)$ we solve the following problem:
%
% $minimize \,\,\,\, m(x_k + p) = f(x_k) + r(x_k)^T J(x_k) p + \frac{1}{2} p^T J(x_k)^T J(x_k) p$
%
% $s.t \,\,\,\, p^T p \le \Delta_k^2$.
%
%% Two search directions
% Since $\frac{\partial}{\partial p} m(x_k + p) = J(x_k)^T r(x_k) + J(x_k)^T J(x_k) p$ and so
% $\left. \frac{\partial}{\partial p} m(x_k + p) \right|_{p = 0} = J(x_k)^T
% r(x_k)$. Therefore we can decrease $m$ from $x_k$ by moving in a
% direction of negative gradient at $p = 0$ with a sufficiently small step 
% i.e, in the direction
%
% $g_k = -J(x_k)^T r(x_k) = -\nabla f(x_k)$
%
% Disregarding the trust region constraint, $m(x_k + p)$ can be minimized
% by choosing $p$ that satisfies:
%
% $\frac{\partial}{\partial p} m(x_k + p) = J(x_k)^T r(x_k) + J(x_k)^T
% J(x_k) p = 0$
%
% Let us denote this solution by $p_{Nk}$ (a Newton-like direction) then:
%
% $J(x_k)^T J(x_k) p_{Nk} = -J(x_k)^T r(x_k)$
%
% Using a singular value decomposition of $J(x_k)$ it is easy to show that
% the solution for $p_{Nk}$ is unique if $J(x_k)$ is non-singular and if
% $J(x_k)$ is singular then there are infinitely many solutions $p_{Nk}$.
%
% Since $\nabla f(x_k) = J(x_k)^T r(x_k)$ we have that $g_k^T \nabla f(x_k)
% \le 0$. So $g_k$ is a descent direction for $f$ at $x_k$.
%
% Now $p_{Nk}^T \nabla f(x_k) = p_{Nk}^T J(x_k)^T r(x_k) = -p_{Nk}^T J(x_k)^T J(x_k)
% p_{Nk} \le 0$. Therefore, $p_{Nk}$ is also a descent direction for $f$ at
% $x_k$.
%
%% 2-D subspace minimization
%
% In order to deal with the trust region constraint, we require that our
% solution $p$ is a linear combination of $g_k$ and $p_{Nk}$ both of which are
% descent directions for $f$ at $x_k$. First we compute $q \times 2$
% orthonormal matrix $Q$ such that:
%
% $[g_k, p_{Nk}] = Q R$
%
% We parameterize the solution as:
%
% $p = Q u$.
%
% The trust region constraint transforms to:
%
% $p^T p = u^T Q^T Q u = u^T u \le \Delta_k^2$. 
%
% The model function transforms to:
%
% $h(u) = m(x_k + Q u) = f(x_k) + r(x_k)^T J(x_k) Q u + \frac{1}{2} u^T Q^T
% J(x_k)^T J(x_k) Q u$. 
%
% Therefore, in terms of $u$ we solve the problem:
%
% $minimize \,\,\,\, h(u) = m(x_k + Q u) = f(x_k) + \bar{g}_k^T u + \frac{1}{2} u^T \bar{B}_k u$
%
% $s.t \,\,\,\, u^T u \le \Delta_k^2$.
%
% where
%
% $\bar{g}_k^T = r(x_k)^T J(x_k) Q$
%
% $\bar{B}_k = Q^T J(x_k)^T J(x_k) Q$
%
% We can solve the 2-D trust region subproblem exactly using
% solveTrustRegionProblemExact.m. Once $u_k$ is found as a solution to this 
% subproblem, $p_k$ is computed using $p_k = Q u_k$. The tentative step is
% then $x_k + p_k = x_k + Q u_k$. The ratio of actual reduction in $f$ to
% predicted reduction using $m$ when moving from $x_k$ to $x_k + p_k$ is:
%
% $\rho_k = \frac{ r(x_k)^T r(x_k) - r(x_k + p_k)^T r(x_k + p_k)}{r(x_k)^T r(x_k) - \{r(x_k) + J(x_k) p\}^T \{r(x_k) + J(x_k) p\}}$
%
% Standard trust region updating strategy is used for updating $\Delta_k$
% and $x_k$ - for example, algorithm 11.5 in Nocedal and Wright.
%
%% Line search Newton method
%
% In the line search method we select a search direction $p_k$ and then
% determine an appropriate step along this direction via backtracking line
% search. To ensure global convergence, we need to ensure that:
%
% $\cos(\theta_k) = -\frac{p_k^T \nabla f(x_k)}{||p_k||_2 ||\nabla
% f(x_k)||_2} \ge \delta$ for some $\delta \in (0,1)$ where
%
% $g_k = -J(x_k)^T r(x_k) = -\nabla f(x_k)$
%
% Whenever the $J(x_k)$ is non-singular, the Newton step is well defined:
%
% $J(x_k)^T J(x_k) p_{Nk} = -J(x_k)^T r(x_k)$
% 
% or
%
% $J(x_k) p_{Nk} = -r(x_k)$
%
% We can take $p_k = p_{Nk}$ and check if $\cos(\theta_k) = -\frac{p_k^T \nabla f(x_k)}{||p_k||_2 ||\nabla
% f(x_k)||_2} \ge \delta$ holds or not. If not, we modify $p_k$ to be equal
% to:
%
% $p_{k} = -(J(x_k)^T J(x_k) + \tau_k I_q )^{-1} \, J(x_k)^T r(x_k)$ where
% $\tau_k$ is chosen so that $\cos(\theta_k) = -\frac{p_k^T \nabla f(x_k)}{||p_k||_2 ||\nabla
% f(x_k)||_2} \ge \delta$ holds. To compute $p_k$ in this case, we can
% exploit the fact that $(J(x_k)^T J(x_k) + \tau_k I_q ) = R^T R$ where:
%
% $[J(x_k); \sqrt{\tau_k} I_q] = Q R$ with $Q$ a $2q \times q$ matrix and
% $R$ a $q \times q$ matrix. 
%
% After selecting $p_k$, a backtracking line search is done as follows:
%
% $\alpha_k = 1.0$
%
% $while \,\, f(x_k + \alpha_k p_k) > f(x_k) + c_1 \alpha_k p_k^T \nabla f(x_k)$
%
% $\alpha_k = \rho \alpha_k$
%
% $end \,\, while$
%
% and the new iterate is:
%
% $x_{k+1} = x_k + \alpha_k p_k$
%
% $c_1 = 10^{-5}$ and $\rho = 0.5$ are sensible choices of $c_1$ and $\rho$.
% $\delta$ can be chosen as $\delta = 10^{-3}$.

        % 1. Check number of input arguments.
        narginchk(2,Inf)

        % 2. Supported methods.
                         Default = 'Default';
                   TrustRegion2D = 'TrustRegion2D';        
        LineSearchModifiedNewton = 'LineSearchModifiedNewton';
                LineSearchNewton = 'LineSearchNewton';
        
        % 3. Create a default options structure.
         dfltTolFun = 1e-6;
           dfltTolX = 1e-12;
        dfltMaxIter = 1000;
        dfltDisplay = 'off';        
        %dfltoptions = statset('TolFun',dfltTolFun ,'TolX'   ,dfltTolX,...
        %                     'MaxIter',dfltMaxIter,'Display',dfltDisplay);
        dfltoptions = struct('TolFun',dfltTolFun ,'TolX'   ,dfltTolX,...
                             'MaxIter',dfltMaxIter,'Display',dfltDisplay);                         
         dfltmethod = TrustRegion2D;        
                         
        % 4. Parse optional args.
          names = {  'Options',   'Method'};
          dflts = {dfltoptions, dfltmethod};
        [options,method] = internal.stats.parseArgs(names,dflts,varargin{:});
        
        % 5. options must be a structure and method must be valid.                 
        assert(isstruct(options));
        internal.stats.getParamVal(method,{Default,TrustRegion2D,LineSearchModifiedNewton,LineSearchNewton},'Method');
               
        % 6. Combine options with dfltOptions to create filled options.
%        options = statset(dfltoptions,options);
        
        % 7. Extract gradTol, stepTol and maxit from the options structure.
        gradTol = options.TolFun;
        stepTol = options.TolX;
          maxit = options.MaxIter;
        
        % 8. Set verbose to true/false.        
        if strcmpi(options.Display,'off')
            verbose = false;
        else
            verbose = true;
        end
        
        % 9. Validate rfun/x0.
        assert(isa(rfun,'function_handle'));
        assert(isnumeric(x0) & isreal(x0) & iscolumn(x0));
        
        % 10. If x0 is 0-by-1, we are done.
        if size(x0,1) == 0
             xstar = x0;
            rxstar = zeros(0,1);
            Jxstar = zeros(0,0);
             cause = 0;
            return;
        end
        
        % 11. Call the appropriate method.
        switch lower(method)            
            case {lower(TrustRegion2D),lower(Default)}                
                [xstar,rxstar,Jxstar,cause] = ...
                               nlesolveTrustRegion2D(rfun,x0,gradTol,stepTol,maxit,verbose);                
            case lower(LineSearchModifiedNewton)                
                [xstar,rxstar,Jxstar,cause] = ...
                    nlesolveLineSearchModifiedNewton(rfun,x0,gradTol,stepTol,maxit,verbose);                
            case lower(LineSearchNewton)                
                [xstar,rxstar,Jxstar,cause] = ...
                            nlesolveLineSearchNewton(rfun,x0,gradTol,stepTol,maxit,verbose);
        end
        
end % end of nlesolve.

function [xstar,rxstar,Jxstar,cause] = nlesolveTrustRegion2D(rfun,x0,gradTol,stepTol,maxit,verbose)

        % 1. Initial r, J, pn and g at x0.
        [rx0,Jx0,pnx0] = rfun(x0);
                   gx0 = -Jx0'*rx0;
        if any(isinf(gx0))
            gx0 = replaceInf(gx0,realmax);
        end
        % TODO validate rx0, Jx0 and pnx0.
        
        % 2. Initial and maximum trust region radius.
        if any(~isfinite(pnx0))
            Delta0 = 1;
        else            
            Delta0 = norm(pnx0);
        end
        DeltaMax = max(Delta0,1e9);          
        
        % 3. Value for eta and tolerance for detecting solutions close to 
        % the boundary of the trust region.
        eta = 5e-4;
        tol = eps(class(x0))^(1/4);
        
        % 4. Convergence flag and iteration counter.        
        found = false;
         iter = 0;
        
        % 5. Initialize x, rx, Jx, pnx, gx and norms for convergence tests.
            x = x0;
           rx = rx0;
           Jx = Jx0;
          pnx = pnx0;
           gx = gx0;
         twonormrx = norm(rx);
         infnormgx = max(abs(gx));
        infnormgx0 = infnormgx;           
           
        % 6. Initialize current trust region radius Delta.   
        Delta = Delta0;        
        
        % 7. Start iterations.
        while( found == false )
            
            % 7.1 Use Newton direction if valid. Otherwise, do a 2D 
            % subspace minimization.
             if ( all(isfinite(pnx)) && norm(pnx) <= Delta )
                 p = pnx;
             else
                % 7.1.1 If Newton direction has NaN/Inf values, do a 1D
                % subspace minimization using gx (i.e., Cauchy point).
                % Otherwise do a 2D subspace minimization using [gx,pnx].
                if any(~isfinite(pnx))
                    [Q,~] = qr(gx,0);
                else
                    [Q,~] = qr([gx,pnx],0);
                end
                
                % 7.1.2 Form gbarx, Bbarx.
                gbarx = -Q'*gx;
                Bbarx = Jx*Q;
                Bbarx = Bbarx'*Bbarx;
                
                % 7.1.3 Solve for step length p.
                if any(isinf(Bbarx(:)))                    
                    Bbarx = replaceInf(Bbarx,realmax);
                end
                if any(isinf(gbarx))                    
                    gbarx = replaceInf(gbarx,realmax);
                end
                u = solveTrustRegionProblemExact(gbarx,Bbarx,Delta);
                p = Q*u;
            end
            twonorms = norm(p);
            
            % 7.2 Compute rho.
            rxp = rfun(x + p);            
            rxtrx = rx'*rx;
            ared = rxtrx - rxp'*rxp;
            pred = rxtrx - norm(rx + Jx*p)^2;
             rho = ared/pred;
            
            % 7.3 Should we take the step?
            if ( rho > eta )
                  stepTaken = true;           
            else
                  stepTaken = false;
            end            
            
            % 7.4 Display convergence info if requested.
            if ( verbose == true )               
                displayConvergenceInfo(iter, twonormrx, infnormgx, twonorms, rho, Delta, stepTaken);
            end
            
            % 7.5 Update current x, rx, Jx, pnx, gx and norms.
            if ( stepTaken == true )
                         x = x + p;
                [rx,Jx,pnx] = rfun(x);
                         gx = -Jx'*rx; 
                  twonormrx = norm(rx);
                  infnormgx = max(abs(gx));
            end
            
            % 7.6 Update trust region radius.
            if ( rho < 0.25 || isnan(rho) )
                Delta = 0.25*Delta;
                % NOTE: The following update can also be used for non
                % linear equations: Delta = 0.25*twonorms.
            else
                if ( rho > 0.75 && abs(twonorms - Delta) <= tol )
                    Delta = min(2*Delta,DeltaMax);
                end
            end   
            
            % 7.7 Convergence tests.
            tau = max( 1, min( twonormrx, infnormgx0 ) );
            if ( infnormgx <= tau * gradTol )
                found = true;
                cause = 0;
            elseif ( twonorms <= stepTol )
                found = true;
                cause = 1;
            elseif ( iter > maxit )
                found = true;
                cause = 2;
            end
            
            % 7.8 Update iteration counter.
            iter = iter + 1;
            
            % 7.9 Display final convergence message.
            if (found == true && verbose == true)
                displayFinalConvergenceMessage(infnormgx, tau, gradTol, twonorms, stepTol, cause);
            end                      
            
        end
        
        % 8. Return requested outputs.
         xstar = x;
        rxstar = rx;
        Jxstar = Jx;

end % end of nlesolveTrustRegion2D.

function [xstar,rxstar,Jxstar,cause] = nlesolveLineSearchModifiedNewton(rfun,x0,gradTol,stepTol,maxit,verbose)

        % 1. Initial r, J, pn and g at x0.
        [rx0,Jx0,pnx0] = rfun(x0);
                   gx0 = -Jx0'*rx0;
                     q = length(x0);
        % TODO validate rx0, Jx0 and pnx0.
        
        % 2. c1, rho and delta.
           c1 = 1e-5;
          rho = 0.5;
        delta = 1e-2;
               
        % 3. Convergence flag and iteration counter.        
        found = false;
         iter = 0;
        
        % 4. Initialize x, rx, Jx, pnx, gx and norms for convergence tests.           
            x = x0;
           rx = rx0;
           Jx = Jx0;
          pnx = pnx0;
           gx = gx0;
         twonormrx = norm(rx);
         infnormgx = max(abs(gx));
        infnormgx0 = infnormgx;
                             
        % 5. Start iterations.
        while( found == false )
            
            % 5.1 Actual value of the cosine of angle between pnx and gx.
            adelta = pnx'*gx/norm(pnx)/norm(gx);
            
            % 5.2 Compute search direction.
            if ( adelta > delta )
                % 5.2.1 Use Newton search direction.
                             p = pnx;
                modifiedNewton = false;
            else
                % 5.2.2 Use Modified Newton search direction.                
                modifiedNewton = true;
                
                done = false;
                beta = norm(Jx,'fro');
                if all(diag(Jx)) > 0
                    tau = 0;
                else
                    tau = beta/2;
                end                
                mag = 2;
                while ( done == false )                    
                    Jxplus = [Jx;sqrt(tau)*eye(q)];
                    rxplus = [rx;zeros(q,1)];
                    p = -(Jxplus \ rxplus);
                    adelta = p'*gx/norm(p)/norm(gx);                    
                    if ( adelta > delta )
                        done = true;
                    else
                        tau = max(tau*mag,beta/2);
                    end                      
                end                
            end
                        
            % 5.3 Backtracking line search.
             alpha = 1.0;
              done = false;
             term1 = 0.5*(rx'*rx);
             term2 = p'*(-gx);
            while( done == false )                
                rxp = rfun(x + alpha*p);                
                isok = 0.5*(rxp'*rxp) - (term1 + c1*alpha*term2) <= 0;                
                if ( isok == true )
                    done = true;
                elseif ( alpha <= stepTol )
                    done = true;
                else
                    alpha = rho*alpha;
                end
            end            
            twonorms = alpha*norm(p);
            
            % 5.4 Display convergence info if requested.
            if ( verbose == true )               
                displayConvergenceInfoLineSearch(iter, twonormrx, infnormgx, twonorms, adelta, alpha, modifiedNewton);
            end
            
            % 5.5 Update current x, rx, Jx, pnx, gx and norms
                      x = x + alpha*p;
            [rx,Jx,pnx] = rfun(x);
                     gx = -Jx'*rx; 
              twonormrx = norm(rx);
              infnormgx = max(abs(gx));                        
            
            % 5.6 Convergence tests.
            tau = max( 1, min( twonormrx, infnormgx0 ) );
            if ( infnormgx <= tau * gradTol )
                found = true;
                cause = 0;
            elseif ( twonorms <= stepTol )
                found = true;
                cause = 1;
            elseif ( iter > maxit )
                found = true;
                cause = 2;
            end
            
            % 5.7 Update iteration counter.
            iter = iter + 1;
            
            % 5.8 Display final convergence message.
            if (found == true && verbose == true)
                displayFinalConvergenceMessage(infnormgx, tau, gradTol, twonorms, stepTol, cause);
            end                      
            
        end
        
        % 6. Return requested outputs.
         xstar = x;
        rxstar = rx;
        Jxstar = Jx;

end % end of nlesolveLineSearchModifiedNewton.

function [xstar,rxstar,Jxstar,cause] = nlesolveLineSearchNewton(rfun,x0,gradTol,stepTol,maxit,verbose)

        % 1. Initial r, J, pn and g at x0.
        [rx0,Jx0,pnx0] = rfun(x0);
                   gx0 = -Jx0'*rx0;
        % TODO validate rx0, Jx0 and pnx0.
        
        % 2. c1, rho.
           c1 = 1e-5;
          rho = 0.5;
               
        % 3. Convergence flag and iteration counter.        
        found = false;
         iter = 0;
        
        % 4. Initialize x, rx, Jx, pnx, gx and norms for convergence tests.           
            x = x0;
           rx = rx0;
           Jx = Jx0;
          pnx = pnx0;
           gx = gx0;
         twonormrx = norm(rx);
         infnormgx = max(abs(gx));
        infnormgx0 = infnormgx;
        modifiedNewton = false;             
        
        % 5. Start iterations.
        while( found == false )
            
            % 5.1 Actual value of the cosine of angle between pnx and gx.
            adelta = pnx'*gx/norm(pnx)/norm(gx);
            
            % 5.2 Always use Newton search direction.
            p = pnx;                        
                        
            % 5.3 Backtracking line search.
             alpha = 1.0;
              done = false;
             term1 = 0.5*(rx'*rx);
             term2 = p'*(-gx);
            while( done == false )                
                rxp = rfun(x + alpha*p);                
                isok = 0.5*(rxp'*rxp) - (term1 + c1*alpha*term2) <= 0;                
                if ( isok == true )
                    done = true;
                elseif ( alpha <= stepTol )
                    done = true;
                else
                    alpha = rho*alpha;
                end
            end            
            twonorms = alpha*norm(p);
            
            % 5.4 Display convergence info if requested.
            if ( verbose == true )               
                displayConvergenceInfoLineSearch(iter, twonormrx, infnormgx, twonorms, adelta, alpha, modifiedNewton);
            end
            
            % 5.5 Update current x, rx, Jx, pnx, gx and norms.
                      x = x + alpha*p;
            [rx,Jx,pnx] = rfun(x);
                     gx = -Jx'*rx; 
              twonormrx = norm(rx);
              infnormgx = max(abs(gx));                        
            
            % 5.6 Convergence tests.
            tau = max( 1, min( twonormrx, infnormgx0 ) );
            if ( infnormgx <= tau * gradTol )
                found = true;
                cause = 0;
            elseif ( twonorms <= stepTol )
                found = true;
                cause = 1;
            elseif ( iter > maxit )
                found = true;
                cause = 2;
            end
            
            % 5.7 Update iteration counter.
            iter = iter + 1;
            
            % 5.8 Display final convergence message.
            if (found == true && verbose == true)
                displayFinalConvergenceMessage(infnormgx, tau, gradTol, twonorms, stepTol, cause);
            end                      
            
        end
        
        % 6. Return requested outputs.
         xstar = x;
        rxstar = rx;
        Jxstar = Jx;

end % end of nlesolveLineSearchNewton.

%=== displayConvergenceInfoLineSearch
function displayConvergenceInfoLineSearch(iter, twonormr, infnormg, twonorms, adelta, alpha, modifiedNewton)
% Helper function to display iteration wise convergence info for line search method.
%
%            iter = iteration number.
%        twonormr = 2-norm of non linear equations.
%        infnormg = infinity norm of the gradient at current iterate.
%        twonorms = two norm of the proposed step size s.
%          adelta = cosine of angle between negative gradient and current search direction.
%           alpha = step size along selected direction.
%  modifiedNewton = true if Newton direction was modified.
%
%  We will display convergence info like this:
%   -----------------------------------------------------------------------------------------
%   ITER     ||r(x)||    ||J(x)^T r(x)||    NORM STEP    cos(theta)     alpha     Mod Newton?
%   -----------------------------------------------------------------------------------------
%      0    +3.533e+01      5.002e+01       4.462e+01    +4.846e-01    1.250e-01      NO
%      1    +7.216e+00      3.721e+01       1.162e+00    +8.681e-01    1.000e+00      NO
%      2    +1.681e+00      1.512e+01       1.600e-01    +9.278e-01    1.000e+00      NO
%      3    +3.322e-02      2.886e-01       3.373e-03    +9.155e-01    1.000e+00      NO
%      4    +1.481e-05      1.286e-04       1.507e-06    +9.150e-01    1.000e+00      NO

    if ( rem(iter,20) == 0 )
        % Display header.
        fprintf('\n');
        fprintf('  -----------------------------------------------------------------------------------------\n');
        fprintf('  ITER     ||r(x)||    ||J(x)^T r(x)||    NORM STEP    cos(theta)     alpha     Mod Newton?\n');
        fprintf('  -----------------------------------------------------------------------------------------\n');
    end
    
    % Display iteration wise convergence info.
    if modifiedNewton == true
        modifiedNewtonString = 'YES';
    else
        modifiedNewtonString = ' NO';
    end
    fprintf('%6d    %+6.3e      %06.3e       %06.3e    %+6.3e    %06.3e     %3s\n', iter, twonormr, infnormg, twonorms, adelta, alpha, modifiedNewtonString);    
end % end of displayConvergenceInfoLineSearch.

%=== displayConvergenceInfo
function displayConvergenceInfo(iter, twonormr, infnormg, twonorms, rho, Delta, stepTaken)
% Helper function to display iteration wise convergence info for trust region method.
%
%         iter = iteration number.
%     twonormr = 2-norm of non linear equations.
%     infnormg = infinity norm of the gradient at current iterate.
%     twonorms = two norm of the proposed step size s.
%          rho = ratio of actual reduction to predicted reduction.
%        Delta = current size of trust region.
%    stepTaken = true if step was taken and otherwise false.
%
%  We will display convergence info like this:
%   -------------------------------------------------------------------------------------
%   ITER     ||r(x)||    ||J(x)^T r(x)||    NORM STEP        RHO       TRUST RAD   ACCEPT
%   -------------------------------------------------------------------------------------
%      0    +4.288e+00      2.945e+01       9.395e-01    +6.328e-01    9.395e-01     YES
%      1    +2.598e+00      2.357e+01       3.635e-01    +9.995e-01    9.395e-01     YES
%      2    +5.808e-02      3.217e-01       3.093e-02    +9.975e-01    9.395e-01     YES
%      3    +2.928e-03      2.624e-02       5.571e-04    +1.000e+00    9.395e-01     YES
%      4    +2.559e-06      2.314e-05       2.467e-07    +1.000e+00    9.395e-01     YES

    if ( rem(iter,20) == 0 )
        % Display header.
        fprintf('\n');
        fprintf('  -------------------------------------------------------------------------------------\n');
        fprintf('  ITER     ||r(x)||    ||J(x)^T r(x)||    NORM STEP        RHO       TRUST RAD   ACCEPT\n');
        fprintf('  -------------------------------------------------------------------------------------\n');
    end
    
    % Display iteration wise convergence info.
    if stepTaken == true
        stepTakenString = 'YES';
    else
        stepTakenString = ' NO';
    end
    fprintf('%6d    %+6.3e      %06.3e       %06.3e    %+6.3e    %06.3e     %3s\n', iter, twonormr, infnormg, twonorms, rho, Delta, stepTakenString);    
end % end of displayConvergenceInfo.

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
    elseif ( cause == 1 )
        % Local minimum possible.
        fprintf([getString(message('stats:classreg:regr:lmeutils:fminqn:Message_LocalMinPossible')),'\n']);
    elseif( cause == 2 )
        % Iteration limit reached.
        fprintf([getString(message('stats:classreg:regr:lmeutils:fminqn:Message_UnableToConverge')),'\n']);
    end
end % end of displayFinalConvergenceMessage.

%=== replaceInf
function B = replaceInf(B,value)
%replaceInf - Replace Inf values in B.
%   B = replaceInf(B,value) takes a matrix or vector B, a scalar value and 
%   replaces +Inf elements in B with abs(value) and -Inf elements in B with
%   -abs(value). The resulting B is returned.

        % (1) Ensure that B is a numeric, matrix and value is a numeric 
        % scalar.
        %assert( isnumeric(B) & ismatrix(B) );
        %assert( isnumeric(value) & isscalar(value) );

        % (2) Get abs(value).
        absvalue = abs(value);
        
        % (3) Find +Inf or -Inf in B.
        isinfB = isinf(B);
        
        % (4) Find +Inf elements in B and replace them with abs(value).
        B(isinfB & B > 0) = absvalue;
        
        % (5) Find -Inf elements in B and replace them with -abs(value).
        B(isinfB & B < 0) = -absvalue;
end % end of replaceInf.

