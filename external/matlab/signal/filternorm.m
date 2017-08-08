function s = filternorm(b,a,pnorm,tol)
%FILTERNORM  Compute the norm of a digital filter.
%   FILTERNORM(B,A) returns the L2-norm of a digital filter given by
%   numerator and denominator coefficients in vectors B and A.  
%
%   FILTERNORM(B,A,PNORM) returns the L2-norm or Linf-norm of a filter.
%   PNORM may be either 2 or inf.  By default, the L2-norm is returned. If
%   computing the L2-norm of an IIR filter, FILTERNORM(...,TOL) will
%   specify the tolerance for greater or less accuracy.  By default,  TOL =
%   1e-8.
%
%   EXAMPLES:
%   % Compute the L2-norm with a tolerance of 1e-10 for an IIR filter
%   [b,a] = butter(5,.5);
%   L2 = filternorm(b,a,2,1e-10); 
%
%   % Compute the infinity norm for a FIR filter
%   b = firpm(30,[.1 .9],[1 1],'Hilbert');
%   Linf = filternorm(b,1,inf);    
%
%   See also ZP2SOS, NORM.

%   Reference:
%     Leland B. Jackson, Digital Filters and Signal Processing, 3rd Ed.
%     Kluwer Academic Publishers, 1996, Chapter 11.

%   Author(s): Dr. Dehner, R. Losada. 
%   Copyright 1988-2004 The MathWorks, Inc.

error(nargchk(2,4,nargin,'struct'));

if nargin == 2, pnorm = 2; end

% Check the input data type. Single precision is not supported.
try
    chkinputdatatype(b,a,pnorm);
catch ME
    throwAsCaller(ME);
end

% Form a df2t filter object
f = dfilt.df2t(b,a);

% Check for stability
if ~isstable(f),
    error(message('signal:filternorm:SignalErrUnstable'));
end

if isinf(pnorm),
    % inf-norm is simply given by the max of the magnitude response
    h = freqz(f,1024);
    s = max(abs(h));
elseif pnorm == 2, 
    % Convert f to transfer function.
    [b,a] = tf(f);
    
    if isfir(f),
        % For a FIR filter, compute 2-norm by simply summing up the 
        % square of the impulse response.        
        s = norm(b,pnorm);
        if nargin == 4,
            warning(message('signal:filternorm:Ignore'));
        end
    else
        % Default tolerance
        if nargin < 4, tol = 1e-8; end
        
        % For a IIR filter, compute 2-norm by approximating the impulse response
        % as finite, alternatively use residues to compute contour integral
        maxradius = max(abs(roots(a)));
        
        % Include an extra check for stability in case numerical roundoff problems
        if maxradius >= 1,
            error(message('signal:filternorm:SignalErrUnstableNumPrec'));
        end
        
        % Determine the number of impulse response points
        N = impzlength(b,a,tol);
        H = impz(b,a,N);
        s = norm(H,pnorm);
    end
end

% [EOF]
