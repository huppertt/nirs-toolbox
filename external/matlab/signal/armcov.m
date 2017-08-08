function [a,e] = armcov( x, p)
%ARMCOV   AR parameter estimation via modified covariance method.
%   A = ARMCOV(X,ORDER) returns the coefficients of the autoregressive (AR)
%   parametric signal model estimate of X using the modified covariance
%   method. The model has order ORDER, and the output array A has ORDER+1
%   columns.  The coefficients along the Nth row of A model the Nth column
%   of X.  If X is a vector then A is a row vector.
%
%   [A,E] = ARMCOV(...) returns the variance estimate E of the white noise
%   input to the AR model.
%
%   % Example:
%   %   Use modified covariance method to estimate the coefficients of an
%   %   autoregressive process given by x(n) = 0.1*x(n-1) -0.8*x(n-2) + 
%   %   w(n).
% 
%   % Generate AR process by filtering white noise
%   a = [1, .1, -0.8];                      % AR coefficients
%   v = 0.4;                                % noise variance
%   w = sqrt(v)*randn(15000,1);             % white noise
%   x = filter(1,a,w);                      % realization of AR process
%   [ar,ec] = armcov(x,numel(a)-1)          % estimate AR model parameters 
% 
%   See also PMCOV, ARCOV, ARBURG, ARYULE, LPC, PRONY.

%   References:
%     [1] S. Lawrence Marple, DIGITAL SPECTRAL ANALYSIS WITH APPLICATIONS,
%              Prentice-Hall, 1987, Chapter 8
%     [2] Steven M. Kay, MODERN SPECTRAL ESTIMATION THEORY & APPLICATION,
%              Prentice-Hall, 1988, Chapter 7

%   Author(s): R. Losada and P. Pacheco
%   Copyright 1988-2014 The MathWorks, Inc.

narginchk(2,2);

[a,e,msg,msgobj] = arparest(x,p,'modified');
if ~isempty(msg), error(msgobj); end

% [EOF] - armcov.m
