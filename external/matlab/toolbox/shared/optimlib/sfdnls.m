function [J,numEvals] = sfdnls(xcurr,valx,Jstr,group,alpha,fun,lb,ub,...
                               finDiffOpts,sizes,finDiffFlags,varargin)
%

%SFDNLS Estimate Jacobian via finite differences
%
%   J = SFDNLS(...) returns an approximation J of the Jacobian matrix of
%   the function 'fun' at the current point xcurr. Dependent on the
%   supplied sparsity structure, Jstr, and coloring information, group,
%   the Jacobian is returned as either sparse or full matrix. 
%
%   In the sparse case:-
%   The vector group indicates how to use sparse finite differencing:
%   group(i) = j means that column i belongs to group (or color) j. Each
%   group (or color) corresponds to a function difference. The input
%   varargin contains the extra parameters (possibly) needed by function
%   'fun'. DiffMinChange and DiffMaxChange indicate, respectively, the
%   minimum and maximum change in variables during the finite difference
%   calculation. A non-empty input alpha overrides the default finite
%   differencing stepsize.
%
%   In the dense case:- 
%   The MATLAB file, finitedifferences, is used to estimate the Jacobian.
%   finitedifferences requires the following inputs to be configured by the
%   caller of this function:-
%   * Lower and upper bounds on the variables (lb, ub).
%   * Options structure (finDiffOpts)
%   * Flags structure (finDiffFlags)
%   * Sizes structure (sizes)
%   More information on the required flags, options and sizes structures 
%   can be found in the MATLAB file help for finitedifferences.
%   
%   [J,numEvals] = SFDNLS(...) returns the number of function evaluations. 

%   Copyright 1990-2011 The MathWorks, Inc.

%
if nargin < 8
   error(message('optimlib:sfdnls:RequiresEightArguments'))
end

scalealpha = false;
x = xcurr(:); % make it a vector
[m,n] = size(Jstr); 
ncol = max(group); 
if isempty(alpha)
    scalealpha = true;
    alpha = repmat(sqrt(eps),ncol,1);
end
J = spones(Jstr);

if ncol < n
   % Estimate Jacobian using sparse finite differences
   for k = 1:ncol
      d = (group == k);
      if scalealpha
         xnrm = norm(x(d));
         xnrm = max(xnrm,1);
         alpha(k) = alpha(k)*xnrm;
      end
      
      % Ensure magnitude of step-size lies within interval 
      % [DiffMinChange, DiffMaxChange]
      alpha(k) = sign(alpha(k))*min(max(abs(alpha(k)),finDiffOpts.DiffMinChange), ...
                                  finDiffOpts.DiffMaxChange);      
      y = x + alpha(k)*d;
      
      xcurr(:) = y;  % reshape for userfunction
      v = feval(fun,xcurr,varargin{:});
      v = v(:);
      
      w = (v-valx)/alpha(k);
      cols = find(d); 
      
      A = sparse(m,n);
      A(:,cols) = J(:,cols);
      J(:,cols) = J(:,cols) - A(:,cols);
      [i,j,val] = find(A);
      [p,ind] = sort(i);
      val(ind) = w(p);
      A = sparse(i,j,full(val),m,n);
      J = J + A;
   end
   numEvals = ncol;
else % ncol == n
   % Estimate Jacobian using dense finite differences
   J = full(J);
   [J,~,~,numEvals] = finitedifferences(x,fun,[],lb,ub,valx,[],[], ...
       1:n,finDiffOpts,sizes,J,[],[],finDiffFlags,[],varargin{:});   
end


