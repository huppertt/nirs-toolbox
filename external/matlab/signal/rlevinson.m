function [R,U,kr,e]=rlevinson(a,efinal)
%RLEVINSON  Reverse Levinson-Durbin Recursion.
%   R=RLEVINSON(A,Efinal) computes the autocorrelation coefficients, R based
%   on the prediction polynomial A and the final prediction error Efinal,
%   using the stepdown algorithm. A should be a minimum phase polynomial
%   and A(1) is assumed to be unity.
%
%   [R,U]=RLEVINSON(...) returns a (P+1)x(P+1) upper triangular matrix, U,
%   that holds the i'th order prediction polynomials Ai, i=1:P, where P
%   is the order of the input polynomial, A.
%
%   The format of this matrix U, is:
%
%         [ 1  a1(1)*  a2(2)* ..... aP(P)  * ]
%         [ 0  1       a2(1)* ..... aP(P-1)* ]
%   U  =  [ .................................]
%         [ 0  0       0      ..... 1        ]
%
%   from which the i'th order prediction polynomial can be extracted
%   using Ai=U(i+1:-1:1,i+1)'. The first row of U contains the
%   conjugates of the reflection coefficients, and the K's may be
%   extracted using, K=conj(U(1,2:end)).
%
%   [R,U,K]=RLEVINSON(...) returns the reflection coefficients in K.
%   [R,U,K,E]=RLEVINSON(...) returns the prediction error of descending
%   orders P,P-1,...,1 in the vector E.
%
%   % Example 1:
%   %   Compute the reflection coefficients of an autoregressive process
%   %   given by x(n) = 0.1*x(n-1) -0.8*x(n-2) + w(n).
%
%   varAR = 0.4;
%   w = sqrt(varAR)*randn(15000,1);
%   a = [1, .1, -0.8];
%   [r,~,k] = rlevinson(a,varAR)    % k is vector of reflection coefficients
%                                   % r is vector of correlation coefficients
%
%   % Compute prediction error values for 1st order, 2nd order predictor
%   p1 = r(1)*(1-k(1)^2)        % 1st order prediction error
%   p2 = p1*(1-abs(k(2))^2)     % 2nd order prediction error, equals varAR
%
%   % Example 2:
%   %   Estimate the spectrum of two sine waves in noise using an AR model.
%   %   Choose the best model order from a group of models returned by the
%   %   reverse Levinson-Durbin recursion.
%
%   Fs = 1000;
%   t = ((0:50e3-1)/Fs).';
%   x = sin(2*pi*50*t) + sin(2*pi*55*t) + 0.2*randn(50e3,1);
%
%   [a,e] = arcov(x,100);               % Estimate AR model parameters
%   [r,u,k] = rlevinson(a,e);
%
%   N = [1,5,25,50,100];
%   nFFT = 8096;
%   P = zeros(nFFT,5);
%
%   % PSD for orders 1,5,25,50,100
%   for idx = 1:numel(N)
%       order = N(idx);
%       ARtest = flipud(u(:,order));
%       P(:,idx) = 1./abs(fft(ARtest,nFFT)).^2;
%   end
%
%   F = (0:1/nFFT:1/2-1/nFFT)*Fs;
%   plot(F, 10*log10(P(1:length(P)/2,:))); grid on
%   legend('Order = 1','Order = 5','Order = 25','Order = 50','Order = 100')
%   xlabel('Frequency Hz'); ylabel('dB')
%   ax = axis; axis([35 70 ax(3:4)])
%
%   See also LEVINSON.

%   References: [1] S. Kay, Modern Spectral Estimation,
%                   Prentice Hall, N.J., 1987, Chapter 6.
%               [2] P. Stoica R. Moses, Introduction to Spectral Analysis
%                   Prentice Hall, N.J., 1997, Chapter 3.
%
%   Author(s): A. Ramasubramanian
%   Copyright 1988-2004 The MathWorks, Inc.

%   Some preliminaries first

% Cast to enforce precision rules 
% This function accepts non-float numeric types but there are no enforced
% rules for the arithmetic in those cases. 
if any([signal.internal.sigcheckfloattype(a,'single','rlevinson','A','allownumeric') ...
    signal.internal.sigcheckfloattype(efinal,'single','rlevinson','Efinal','allownumeric')])
    a = single(a);
    efinal = single(efinal);
end

a = a(:);                   % Convert to a column vector if not already so
if a(1)~=1,
    warning(message('signal:rlevinson:InvalidParam'));
    a = a/a(1);
end

p = length(a);

if p < 2
    error(message('signal:rlevinson:InvalidDimensions'));
end

% This matrix will have the prediction polynomials % of orders 1:p
% Cast to enforce precision rules
if isa(a,'single')
  U = zeros(p,p,'single');
else
  U = zeros(p,p);
end

U(:,p) = conj(a(end:-1:1)); % Prediction coefficients of order p

p = p-1;

% First we find the prediction coefficients of smaller orders and form the
% Matrix U

% Initialize the step down
e(p) = efinal;             % Prediction error of order p

% Step down
for k=p-1:-1:1
    [a,e(k)] = levdown(a,e(k+1));
    U(:,k+1) = [conj(a(end:-1:1).');zeros(p-k,1)];
end

e0 = e(1)/(1-abs(a(2)^2)); % Because a[1]=1 (true polynomial)
U(1,1) = 1;                % Prediction coefficient of zeroth order
kr = conj(U(1,2:end));     % The reflection coefficients
kr = kr.';                 % To make it into a column vector

% Once we have the matrix U and the prediction error at various orders, we can
% use this information to find the autocorrelation coefficients.

% Initialize recursion
R0 = e0;                   % To take care of the zero indexing problem
R(1) = -conj(U(1,2))*R0;   % R[1]=-a1[1]*R[0]

% Actual recursion
for k = 2:1:p
    R(k) = -sum(conj(U(k-1:-1:1,k)).*R(end:-1:1).')-kr(k)*e(k-1); %#ok<AGROW>
end

% Include R(0) and make it a column vector. Note the dot transpose

R = [R0 R].';

% [EOF] rlevinson.m
