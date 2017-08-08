function varargout = powerest(~,~,~) %#ok<STOUT>
%POWEREST   Computes the powers and frequencies of sinusoids.
%   
%   POWEREST is not recommended.  Use <a href="matlab:help rootmusic">rootmusic</a> and <a href="matlab:help rooteig">rooteig</a> instead.
%
%   POW = POWEREST(H,X) returns the vector POW containing the estimates
%   of the powers of the complex sinusoids contained in the data
%   represented by X.  H must be a <a href="matlab:help spectrum.music">music</a> or <a href="matlab:help spectrum.eigenvector">eigenvector</a> estimator. 
%
%   X can be a vector or a matrix. If it's a vector it is a signal, if it's
%   a matrix it may be either a data matrix such that X'*X=R, or a
%   correlation matrix R.  How X is interpreted depends on the spectral
%   estimator's (H) input type, which can be any one of the following:
%       'Vector'  (default)
%       'DataMatrix'
%       'CorrelationMatrix'
%
%   [POW,W] = POWEREST(...) returns in addition a vector of frequencies W
%   of the sinusoids contained in X.  W is in units of rad/sample.
%
%   [POW,F] = POWEREST(...,Fs) uses the sampling frequency Fs in the
%   computation and returns the vector of frequencies, F, in Hz.
%
%   EXAMPLE:
%      n = 0:99;   
%      s = exp(1i*pi/2*n)+2*exp(1i*pi/4*n)+exp(1i*pi/3*n)+randn(1,100);  
%      H = spectrum.music(3);
%      [P,W] = powerest(H,s);

%   Author(s): P. Pacheco
%   Copyright 1988-2012 The MathWorks, Inc.

% Help for powerest.m

% [EOF]
