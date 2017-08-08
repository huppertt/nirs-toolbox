function [xf,fval,exitflag,output] = lclfminbnd(funfcn,ax,bx,options,varargin)
%LCLFMINBND Scalar bounded nonlinear function minimization.  
%
%   NOTE: This is a copy of the toolbox\matlab\funfun\fminbnd.m which
%   shipped in R12.1.  There were changes made to the R13 version of this 
%   file which affected the CHEB1ORD, CHEB2ORD, BUTTORD, and ELLIPORD functions 
%   in the Signal Processing Toolbox (g142835).
%
%   X = FMINBND(FUN,x1,x2) starts at X0 and finds a local minimizer X of the
%   function FUN in the interval x1 < X < x2. FUN accepts scalar input X and returns 
%   a scalar function value F evaluated at X.  
%
%   X = FMINBND(FUN,x1,x2,OPTIONS) minimizes with the default optimization
%   parameters replaced by values in the structure OPTIONS, created with
%   the OPTIMSET function.  See OPTIMSET for details.  FMINBND uses these options: 
%   Display, TolX, MaxFunEval and MaxIter. 
%
%   X = FMINBND(FUN,x1,x2,OPTIONS,P1,P2,...) provides for additional
%   arguments, which are passed to the objective function, FUN(X,P1,P2,...).
%   (Use OPTIONS = [] as a place holder if no options are set.)
%
%   [X,FVAL] = FMINBND(...) also returns the value of the objective function,
%   FVAL, computed in FUN, at X.
%
%   [X,FVAL,EXITFLAG] = FMINBND(...) also returns a string EXITFLAG that 
%   describes the exit condition of FMINBND.  
%   If EXITFLAG is:
%     1 then FMINBND converged with a solution X based on OPTIONS.TolX
%     0 then the maximum number of function evaluations was reached
%   
%   [X,FVAL,EXITFLAG,OUTPUT] = FMINBND(...) also returns a structure
%   OUTPUT with the number of iterations taken in OUTPUT.iterations.
%
%   Examples
%     FUN can be specified using @:
%        X = fminbnd(@cos,3,4)
%      computes pi to a few decimal places and gives a message upon termination.
%        [X,FVAL,EXITFLAG] = fminbnd(@cos,3,4,optimset('TolX',1e-12,'Display','off')) 
%      computes pi to about 12 decimal places, suppresses output, returns the  
%      function value at x, and returns an EXITFLAG of 1.
%
%     FUN can also be an anonymous function handle:
%        f = @(x) sin(x)+3;
%        x = fminbnd(f,2,5);
%
%   See also OPTIMSET, FMINSEARCH, FZERO, @.

%   Reference: "Computer Methods for Mathematical Computations",
%   Forsythe, Malcolm, and Moler, Prentice-Hall, 1976.

%   Original coding by Duane Hanselman, University of Maine.
%   Copyright 1984-2009 The MathWorks, Inc.


defaultopt = struct('Display','notify',...
     'MaxFunEvals',500,'MaxIter',500,'TolX',1e-4);
  
% If just 'defaults' passed in, return the default options in X
if nargin==1 & nargout <= 1 & isequal(funfcn,'defaults')
   xf = defaultopt;
   return
end

if nargin < 3
   error(message('signal:lclfminbnd:Nargchk'))
end

% initialization
if nargin<4, 
   options = []; 
end

printtype = optimget(options,'Display',defaultopt,'fast');
tol = optimget(options,'TolX',defaultopt,'fast');
maxfun = optimget(options,'MaxFunEvals',defaultopt,'fast');
maxiter = optimget(options,'MaxIter',defaultopt,'fast');

% fun evals are the same as iterations, so just check the minimum
maxcount = min(maxfun,maxiter);

switch printtype
case 'notify'
   print = 1;
case {'none','off'}
   print = 0;
case 'iter'
   print = 3;
case 'final'
   print = 2;
otherwise
   print = 1;
end
% checkbounds
if ax > bx
   exitflag = -1;
   output = []; xf=[]; fval = [];
   if print > 0
        msg=sprintf([getString(message('signal:lclfminbnd:ExitingDueToInfeasibility'))]);
        disp(msg)
   end
   return
end

% Assume we'll converge
exitflag = 1;

header = ' Func-count     x          f(x)         Procedure';
step='       initial';

% Convert to inline function as needed.
funfcn = fcnchk(funfcn,length(varargin));

seps = sqrt(eps);
c = 0.5*(3.0 - sqrt(5.0));
a = ax; b = bx;
v = a + c*(b-a);
w = v; xf = v;
d = 0.0; e = 0.0;
x= xf; fx = feval(funfcn,x,varargin{:}); 
num = 1;
fmin_data = [ 1 xf fx ];
if print > 2
   disp(' ')
   disp(header)
   disp([sprintf('%5.0f   %12.6g %12.6g ',fmin_data), step]) 
end
fv = fx; fw = fx;
xm = 0.5*(a+b);
tol1 = seps*abs(xf) + tol/3.0;   
tol2 = 2.0*tol1;

% Main loop
while ( abs(xf-xm) > (tol2 - 0.5*(b-a)) ) 
   gs = 1;
   % Is a parabolic fit possible
   if abs(e) > tol1
      % Yes, so fit parabola
      gs = 0;
      r = (xf-w)*(fx-fv);
      q = (xf-v)*(fx-fw);
      p = (xf-v)*q-(xf-w)*r;
      q = 2.0*(q-r);
      if q > 0.0,  p = -p; end
      q = abs(q);
      r = e;  e = d;
      
      % Is the parabola acceptable
      if ( (abs(p)<abs(0.5*q*r)) & (p>q*(a-xf)) & (p<q*(b-xf)) )
         
         % Yes, parabolic interpolation step
         d = p/q;
         x = xf+d;
         step = '       parabolic';
         
         % f must not be evaluated too close to ax or bx
         if ((x-a) < tol2) | ((b-x) < tol2)
            si = sign(xm-xf) + ((xm-xf) == 0);
            d = tol1*si;
         end
      else
         % Not acceptable, must do a golden section step
         gs=1;
      end
   end
   if gs
      % A golden-section step is required
      if xf >= xm, e = a-xf;    else, e = b-xf;  end
      d = c*e;
      step = '       golden';
   end
   
   % The function must not be evaluated too close to xf
   si = sign(d) + (d == 0);
   x = xf + si * max( abs(d), tol1 );
   fu = feval(funfcn,x,varargin{:});  
   num = num+1;
   fmin_data = [num x fu];
   if print > 2
      disp([sprintf('%5.0f   %12.6g %12.6g ',fmin_data), step]) 
   end
   
   % Update a, b, v, w, x, xm, tol1, tol2
   if fu <= fx
      if x >= xf, a = xf; else, b = xf; end
      v = w; fv = fw;
      w = xf; fw = fx;
      xf = x; fx = fu;
   else % fu > fx
      if x < xf, a = x; else,b = x; end
      if ( (fu <= fw) | (w == xf) )
         v = w; fv = fw;
         w = x; fw = fu;
      elseif ( (fu <= fv) | (v == xf) | (v == w) )
         v = x; fv = fu;
      end
   end
   xm = 0.5*(a+b);
   tol1 = seps*abs(xf) + tol/3.0; tol2 = 2.0*tol1;
   if num >= maxcount
       exitflag = 0;
       output.iterations=num;
       output.funcCount=num;
       fval = fx;
       if print > 0     
           terminate(x,exitflag,fval,maxfun,maxiter,tol,print);
       end
       return
   end
end
fval = fx;
output.iterations = num;
output.funcCount = num;
output.algorithm = 'golden section search, parabolic interpolation';
if print > 0
   terminate(x,exitflag,fval,maxfun,maxiter,tol,print);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function terminate(x,exitflag,finalf,maxfun,maxiter,tol,print)

switch exitflag
case 1
    if print > 1 % only print success if not 'off' or 'notify'
        convmsg1 = getString(message('signal:lclfminbnd:OptimizationTerminatedSuccessfully', sprintf('%e',tol)));
        disp(convmsg1)
    end
case 0
   if maxfun <= maxiter
      disp(' ')
      warning(message('signal:lclfminbnd:DidNotConvergeFunctionEvals', 'MaxFunEvals'));
   else
      disp(' ')
      warning(message('signal:lclfminbnd:DidNotConvergeIters', 'MaxIter'));
      msg = sprintf(['         ' getString(message('signal:lclfminbnd:CurrentFunctionValue')) ': %f \n'], finalf);
      disp(msg)
   end
   
case -1
   ;
otherwise
   ;
end
