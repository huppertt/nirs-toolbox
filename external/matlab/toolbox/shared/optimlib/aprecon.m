function[RPCMTX,ppvec] = aprecon(A,upperbandw,DM,DG,varargin)
%

%APRECON Banded preconditioner function for least-squares problems.
%
%   [RPCMTX,PPVEC] = APRECON(A,UPPERBW,DM,DG) produces the sparse
%   nonsingular upper triangular matrix RPCMTX such that
%   RPCMTX'*RPCMTX is a preconditioner for the
%   matrix M = DM*(A'*A)*DM + DG, where DM is a positive
%   diagonal matrix, DG is a non-negative diagonal matrix,
%   and A is sparse rectangular matrix with more rows than columns.
%   PPVEC is the associated permutation (row) vector and UPPERBW
%   specifies the upperbandwidth of RPCMTX.

%   Default preconditioner for SLLSBOX and SNLS.

%   Copyright 1990-2011 The MathWorks, Inc.

% Initialization
[rows,cols] = size(A);
n = length(DM);
RPCMTX = speye(n); 
ppvec = (1:n);
% In case "A" isn't really A, but something else to use with JacobMult function.
if ~isnumeric(A) || isempty(A) || (rows < cols) || ~isequal(n,cols)
    % A is not the right size; ignore requested bandwidth and compute
    % diagonal preconditioner based only on DM and DG.
    d1 = full(diag(DM)); 
    d2 = full(diag(DG)); 
    dd = sqrt(d1.*d1 + abs(d2));
    RPCMTX = sparse(1:n,1:n,dd);
    return
end

epsi = sqrt(eps);
if nargin < 2 || isempty(upperbandw)
   upperbandw = 0; 
end

% Form matrix M
TM = A*DM;

% Determine factor of preconditioner.
if upperbandw == 0 % Diag. preconditioner based on column 2-norms
   M = TM'*TM + DG;
   % M cannot be a row vector so no problem with sum.
   dnrms = sqrt(sum(M.*M))';
   d = max(sqrt(dnrms),epsi);
   RPCMTX = sparse(1:n,1:n,full(d));
   ppvec = (1:n);
   % But what if RPCMTX is singular?
elseif upperbandw >= n-1 % Attempt sparse QR
   dgg = sqrt(diag(DG)); 
   dgg = full(dgg);
   DDG = sparse(1:n,1:n,dgg);
   TM = [TM;DDG];
   p = colamd(TM);
   RPCMTX = qr(TM(:,p)); 
   RPCMTX = RPCMTX(1:n,1:n);
   ppvec = p;
   
   %    Modify for singularity?
   mdiag = min(abs(diag(RPCMTX)));
   lambda = 1;
   while mdiag < sqrt(eps);
      TM = [A*DM; DDG + lambda*speye(n)];
      lambda = 4*lambda;
      p = colamd(TM);
      RPCMTX = qr(TM(:,p));
      RPCMTX = RPCMTX(1:n,1:n);
      ppvec = p;
      mdiag = min(abs(diag(RPCMTX)));
   end
   
elseif (upperbandw > 0) && (upperbandw < n-1) % Band approximation.
   M = TM'*TM + DG;
   p = (1:n);
   M = tril(triu(M(p,p),-upperbandw),upperbandw);
   RPCMTX = sparse(n,n);
   [RPCMTX,info] = chol(M);
   lambda = 1;
   while info > 0
      M = M + lambda*speye(n);
      RPCMTX = sparse(n,n);
      [RPCMTX,info] = chol(M);
      lambda = 4*lambda;
   end
   ppvec = p;
else
   error(message('optimlib:aprecon:InvalidUpperbandw'));
end


