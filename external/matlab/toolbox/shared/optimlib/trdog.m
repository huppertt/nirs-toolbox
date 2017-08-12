function[s,snod,qpval,posdef,pcgit,Z] = trdog(x,g,H,D,delta,dv,...
   mtxmpy,pcmtx,pcoptions,tol,kmax,theta,l,u,Z,dnewt,preconflag,varargin)
%

%TRDOG Reflected (2-D) trust region trial step (box constraints)
%
% [s,snod,qpval,posdef,pcgit,Z] = TRDOG(x,g,H,D,delta,dv,...
%                 mtxmpy,pcmtx,pcoptions,tol,theta,l,u,Z,dnewt,preconflag);
%
%   Determine the trial step `s', an approx. trust region solution.
%   `s' is chosen as the best of 3 steps: the scaled gradient
%   (truncated to  maintain strict feasibility),
%   a 2-D trust region solution (truncated to remain strictly feas.),
%   and the reflection of the 2-D trust region solution,
%   (truncated to remain strictly feasible).
%
%   The 2-D subspace (defining the trust region problem) is defined 
%   by the scaled gradient direction and a CG process (returning
%   either an approximate Newton step of a direction of negative curvature.
%   Driver functions are: SNLS, SFMINBX
%   SNLS actually calls TRDOG with the Jacobian matrix (and a special 
%   Jacobian-matrix multiply function in MTXMPY).

%   Copyright 1990-2011 The MathWorks, Inc.

% Initialization
n = length(g);  
pcgit = 0; 
grad = D*g;
DM = D; 
DG = sparse(1:n,1:n,full(abs(g).*dv));
posdef = 1; 
pcgit = 0; 
tol2 = sqrt(eps);
v1 = dnewt; 
qpval1 = inf; 
qpval2 = inf; 
qpval3 = inf;

% DETERMINE A 2-DIMENSIONAL SUBSPACE
if isempty(Z)
   if isempty(v1)
      switch preconflag
      case 'hessprecon'
         % preconditioner based on H, no matter what it is
         [R,permR] = feval(pcmtx,H,pcoptions,DM,DG,varargin{:});
      case 'jacobprecon'
         [R,permR] = feval(pcmtx,H,pcoptions,DM,DG,varargin{:});
      otherwise
         error(message('optimlib:trdog:InvalidPreconflag'));
      end
      % We now pass kmax in from calling function
      %kmax = max(1,floor(n/2));
      if tol <= 0, 
         tol = .1; 
      end 
   
      [v1,posdef,pcgit] = pcgr(DM,DG,grad,kmax,tol,...
         mtxmpy,H,R,permR,preconflag,pcoptions,varargin{:});
   end
   if norm(v1) > 0
      v1 = v1/norm(v1);
   end
   Z(:,1) = v1;
   if n > 1
      if (posdef < 1)
         v2 = D*sign(grad); 
         if norm(v2) > 0
            v2 = v2/norm(v2);
         end
         v2 = v2 - v1*(v1'*v2); 
         nrmv2 = norm(v2);
         if nrmv2 > tol2
            v2 = v2/nrmv2; 
            Z(:,2) = v2;
         end
      else
         if norm(grad) > 0
            v2 = grad/norm(grad);
         else
            v2 = grad;
         end
         v2 = v2 - v1*(v1'*v2); 
         nrmv2 = norm(v2);
         if nrmv2 > tol2
            v2 = v2/nrmv2; 
            Z(:,2) = v2; 
         end
      end
   end
end

%  REDUCE TO THE CHOSEN SUBSPACE
W = DM*Z;  
switch preconflag
case 'hessprecon'
   WW = feval(mtxmpy,H,W,varargin{:}); 
case 'jacobprecon'
   WW = feval(mtxmpy,H,W,0,varargin{:}); 
otherwise
   error(message('optimlib:trdog:InvalidPreconflag'));
end

W = DM*WW;
MM = full(Z'*W + Z'*DG*Z); 
rhs = full(Z'*grad);

%  Determine 2-D TR soln
[st,qpval,po,fcnt,lambda] = trust(rhs,MM,delta);
ss = Z*st;  
s = abs(diag(D)).*ss; 
s = full(s); 
ssave = s;
sssave = ss; 
stsave = st;

% Check direction for NaNs
if any(isnan(s))
   error(message('optimlib:trdog:NaNInStep'))
end

% Truncate the TR solution?
arg = (abs(s) > 0);
% No truncation if s is an all zero vector
if ~any(arg)   
   alpha = 1;
   mmdis = 1;
   ipt = [];  % Set to empty if all-zero step so that it won't be undefined
else
   dis = max((u(arg)-x(arg))./s(arg), (l(arg)-x(arg))./s(arg));
   [mmdis,ipt] = min(dis);  
   mdis = theta*mmdis;
   alpha = min(1,mdis);
end
s = alpha*s; 
st = alpha*st; 
ss = full(alpha*ss);
qpval1 = rhs'*st + (.5*st)'*MM*st;
if n > 1 
   %   Evaluate along the reflected direction?
   qpval3 = inf; 
   ssssave = mmdis*sssave;
   if norm(ssssave) < .9*delta
      r = mmdis*ssave; 
      nx = x+r;
      stsave = mmdis*stsave;
      qpval0 = rhs'*stsave + (.5*stsave)'*MM*stsave;
      switch preconflag
      case 'hessprecon'
         ng = feval(mtxmpy,H,r,varargin{:}); 
      case 'jacobprecon'
         ng = feval(mtxmpy,H,r,0,varargin{:}); 
      otherwise
         error(message('optimlib:trdog:InvalidPreconflag'));
      end
      
      ng = ng + g; 
      ngrad = D*ng;
      ngrad = ngrad + DG*ssssave;
      
      %      nss is the reflected direction
      nss = sssave; 
      nss(ipt) = -nss(ipt); 
      ZZ(:,1) = nss/norm(nss);
      W = DM*ZZ; 
      
      switch preconflag
      case 'hessprecon'
         WW = feval(mtxmpy,H,W,varargin{:}); 
      case 'jacobprecon'
         WW = feval(mtxmpy,H,W,0,varargin{:}); 
      otherwise
         error(message('optimlib:trdog:InvalidPreconflag'));
      end
      
      
      W = DM*WW;
      MM = full(ZZ'*W + ZZ'*DG*ZZ);
      nrhs=full(ZZ'*ngrad);
      [nss,tau] = quad1d(nss,ssssave,delta); 
      nst = tau/norm(nss);
      ns = abs(diag(D)).*nss; 
      ns = full(ns);

      % Check direction for NaNs
      if any(isnan(ns))
         error(message('optimlib:trdog:NaNInRflctdStep'))
      end
      
      % Truncate the reflected direction?
      arg = (abs(ns) > 0);
      % No truncation if s is zero length
      if ~any(arg) 
         alpha = 1;
      else
         dis = max((u(arg)-nx(arg))./ns(arg), (l(arg)-nx(arg))./ns(arg));
         mdis = min(dis);  
         mdis = theta*mdis;
         alpha = min(1,mdis); 
      end
      ns = alpha*ns; 
      nst = alpha*nst; 
      nss = full(alpha*nss);
      qpval3 = qpval0 +  nrhs'*nst + (.5*nst)'*MM*nst;
   end
   
   %   Evaluate along gradient direction
   gnorm = norm(grad);
   ZZ(:,1) = grad/(gnorm + (gnorm == 0)); % Protect against norm of 0
   W = DM*ZZ; 
   
   switch preconflag
   case 'hessprecon'
      WW = feval(mtxmpy,H,W,varargin{:}); 
   case 'jacobprecon'
      WW = feval(mtxmpy,H,W,0,varargin{:}); 
   otherwise
      error(message('optimlib:trdog:InvalidPreconflag'));
   end
   
   
   W = DM*WW;
   MM = full(ZZ'*W + ZZ'*DG*ZZ); 
   rhs = full(ZZ'*grad);
   [st,qpval,po,fcnt,lambda] = trust(rhs,MM,delta);
   ssg = ZZ*st; 
   sg = abs(diag(D)).*ssg; 
   sg = full(sg);
   
   % Check direction for NaNs
   if any(isnan(sg))
      error(message('optimlib:trdog:NaNInGrad'))
   end
   
   %   Truncate the gradient direction?
   arg = (abs(sg) > 0);
   if ~any(arg)
       % No truncation if s is zero length
      alpha = 1;
   else
      dis = max((u(arg)-x(arg))./sg(arg), (l(arg)-x(arg))./sg(arg));
      mdis = min(dis); 
      mdis = theta*mdis;
      alpha = min(1,mdis); 
   end
   sg = alpha*sg; 
   st = alpha*st; 
   ssg = full(alpha*ssg);
   qpval2 = rhs'*st + (.5*st)'*MM*st;
end

% Choose the best of s, sg, ns.
if qpval2 <= min(qpval1,qpval3)
   qpval = qpval2; 
   s = sg; 
   snod = ssg;
elseif qpval1 <= min(qpval2,qpval3)
   qpval = qpval1; 
   snod = ss;
else
   qpval = qpval3; 
   s = ns + r; 
   snod = nss + ssssave;
end

%-----------------------------------------------------------
function[nx,tau] = quad1d(x,ss,delta)
%QUAD1D	1D quadratic zero finder for trust region step
%
% [nx,tau] = quad1d(x,ss,delta) tau is min(1,step-to-zero)
% of a 1-D quadratic ay^2 + b*y + c.
% a = x'*x; b = 2*(ss'*x); c = ss'*ss-delta^2). nx is the
% new x value, nx = tau*x;

% Algorithm:
% numer = -(b + sign(b)*sqrt(b^2-4*a*c));
% root1 = numer/(2*a);
% root2 = c/(a*root1);   % because root2*root1 = (c/a);

a = x'*x;
b = 2*(ss'*x); 
c = ss'*ss-delta^2;

numer = -(b + sign(b)*sqrt(b^2-4*a*c));
r1 = numer/(2*a);
r2 = c/(a*r1);

tau = max(r1,r2); 
tau = min(1,tau);
if tau <= 0, 
   error(message('optimlib:trdog:SqrRt')); 
end
nx = tau*x;



