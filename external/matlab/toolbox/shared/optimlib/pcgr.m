function[p,posdef,k] = pcgr(DM,DG,g,kmax,tol,mtxmpy,H,R,pR,callerflag,pcoptions,varargin)
%

%PCGR	Preconditioned conjugate gradients
%
% [p,posdef,k] = PCGR(DM,DG,g,kmax,tol,mtxmpy,H,R,pR) apply
% a preconditioned conjugate gradient procedure to the quadratic
% 
%         q(p) = .5p'Mp + g'p, where
%
% M = DM*H*DM + DG. kmax is a bound on the number of permitted
% CG-iterations, tol is a stopping tolerance on the residual (default
% is tol = .1), mtxmpy is the function that computes products
% with the Hessian matrix H,  
% and R is the cholesky factor of the preconditioner (transpose) of
% M. So, R'R approximates M(pR,pR), where pR is a permutation vector.
% On output p is the computed direction, posdef = 1 implies
% only positive curvature (in M) has been detected; posdef = 0
% implies p is a direction of negative curvature (for M).
% Output parameter k is the number of CG-iterations used (which
% corresponds to the number of multiplications with H).
%

%   Copyright 1990-2011 The MathWorks, Inc.

% Initializations.
n = length(DG);
r = -g; 
p = zeros(n,1); 

% Precondition .
z = preproj(r,R,pR);
znrm = norm(z); 
stoptol = tol*znrm;
inner2 = 0; 
inner1 = r'*z;
posdef = 1;

kmax = max(kmax,1);  % kmax must be at least 1
% PRIMARY LOOP.
for k = 1:kmax
   if k==1
      d = z;
   else
      beta = inner1/inner2;
      d = z + beta*d;
   end
   ww = DM*d;
   switch callerflag
   case 'hessprecon'
      w = feval(mtxmpy,H,ww,varargin{:});
   case 'jacobprecon'
      w = feval(mtxmpy,H,ww,0,varargin{:});
   otherwise
      error(message('optimlib:pcgr:InvalidCallingFcn')) 
   end
   ww = DM*w +DG*d;
   denom = d'*ww;
   if denom <= 0
      if norm(d) == 0
         p = d;
      else      
         p = d/norm(d);
      end
      posdef = 0;
      break
   else
      alpha = inner1/denom;
      p = p + alpha*d;
      r = r - alpha*ww;
   end
   z = preproj(r,R,pR);
   
   % Exit?
   if norm(z) <= stoptol 
      break; 
   end
   inner2 = inner1;
   inner1 = r'*z;
end

% If the Newton step is computed via the direct factorization, then it is
% "exact" and no actual PCG iterations took place. This is done for display
% purposes only.
if isinf(pcoptions)
    k = 0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function w = preproj(r,RPCMTX,ppvec)
%PREPROJ Apply preconditioner
%
% w = preproj(r,RPCMTX,ppvec) Apply a preconditioner to vector r.
% The conceptual preconditioner is H(ppvec,ppvec) = RPCMTX'*RPCMTX

% Initialization
n = length(r);

if nargin < 3 || isempty(ppvec)
   ppvec = (1:n); 
   if nargin < 2 || isempty(RPCMTX) 
      RPCMTX = speye(n); 
   end
end 

%  Precondition
wbar = RPCMTX'\r(ppvec);
w(ppvec,1) = RPCMTX\wbar;




