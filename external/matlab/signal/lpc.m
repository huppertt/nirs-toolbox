function [a,e]=lpc(x,N)
%LPC  Linear Predictor Coefficients.
%   A = LPC(X,N) finds the coefficients, A=[ 1 A(2) ... A(N+1) ], of an Nth
%   order forward linear predictor.
%
%      Xp(n) = -A(2)*X(n-1) - A(3)*X(n-2) - ... - A(N+1)*X(n-N)
%
%   such that the sum of the squares of the errors
%
%      err(n) = X(n) - Xp(n)
%
%   is minimized.  X can be a vector or a matrix.  If X is a matrix
%   containing a separate signal in each column, LPC returns a model
%   estimate for each column in the rows of A.  N specifies the order of
%   the polynomial A(z) which must be a positive integer.  N must be less
%   or equal to the length of X.  If X is a matrix, N must be less or equal
%   to the length of each column of X.
%
%   If you do not specify a value for N, LPC uses a default N =
%   length(X)-1.
%
%   [A,E] = LPC(X,N) returns the variance (power) of the prediction error.
%
%   LPC uses the Levinson-Durbin recursion to solve the normal equations
%   that arise from the least-squares formulation.  This computation of the
%   linear prediction coefficients is often referred to as the
%   autocorrelation method.
%
%   % Example:
%   %   Estimate a data series using a third-order forward predictor, and 
%   %   compare to the original signal.
%   
%   %   Create signal data as the output of an autoregressive process
%   %   driven by white noise. Use the last 4096 samples of the AR process
%   %   output to avoid start-up transients:
%   randn('state',0); 
%   noise = randn(50000,1);  % Normalized white Gaussian noise
%   x = filter(1,[1 1/2 1/3 1/4],noise);
%   x = x(45904:50000);
%   a = lpc(x,3);
%   est_x = filter([0 -a(2:end)],1,x);  % Estimated signal
%   e = x - est_x;                      % Prediction error
%   [acs,lags] = xcorr(e,'coeff');      % ACS of prediction error
%   
%   %   Compare the predicted signal to the original signal
%   plot(1:97,x(4001:4097),1:97,est_x(4001:4097),'--');
%   title('Original Signal vs. LPC Estimate');
%   xlabel('Sample Number'); ylabel('Amplitude'); grid;
%   legend('Original Signal','LPC Estimate')
%
%   %   Look at the autocorrelation of the prediction error.
%   figure; plot(lags,acs); 
%   title('Autocorrelation of the Prediction Error');
%   xlabel('Lags'); ylabel('Normalized Value'); grid;
%
%   See also LEVINSON, ARYULE, PRONY, STMCB.

%   Author(s): T. Krauss, 9-21-93
%   Modified:  T. Bryan 11-14-97
%   Copyright 1988-2004 The MathWorks, Inc.

error(nargchk(1,2,nargin,'struct'))

if isempty(x)
    error(message('signal:lpc:Empty'));
end

[m,n] = size(x);
if (n>1) && (m==1)
	x = x(:);
	[m,n] = size(x);
end

if nargin < 2,
    N = m-1; 
elseif N < 0,
    % Check for N positive
    error(message('signal:lpc:negativeOrder'));
end

% Check the input data type. Single precision is not supported.
try
    chkinputdatatype(x,N);
catch ME
    throwAsCaller(ME);
end

if (N > m),
    error(message('signal:lpc:orderTooLarge', 'X must be a vector with length greater or equal to the prediction order.', 'If X is a matrix, the length of each column must be greater or equal to', 'the prediction order.'));
end

% Compute autocorrelation vector or matrix
X = fft(x,2^nextpow2(2*size(x,1)-1));
R = ifft(abs(X).^2);
R = R./m; % Biased autocorrelation estimate

[a,e] = levinson(R,N);

% Return only real coefficients for the predictor if the input is real
for k = 1:n,
    if isreal(x(:,k))
        a(k,:) = real(a(k,:));
    end
end
