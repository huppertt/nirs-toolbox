function varargout = rooteig(varargin)
%ROOTEIG   Computes the frequencies and powers of sinusoids via the
%          eigenvector algorithm.
%   W = ROOTEIG(X,P) returns the vector of frequencies W of the complex
%   sinusoids contained in signal vector X.  W is in units of rad/sample.
%   P is the number of complex sinusoids in X. If X is a data matrix,
%   each row is interpreted as a separate sensor measurement or trial.
%   In this case, X must have a number of columns larger than P.  You can
%   use the function CORRMTX to generate data matrices to be used here.  
%
%   W = ROOTEIG(R,P,'corr') returns the vector of frequencies W, for a
%   signal whose correlation matrix estimate is given by the positive
%   definite matrix R. Exact conjugate-symmetry of R is ensured by forming
%   (R+R')/2 inside the function. The number of rows or columns of R must
%   be greater than P.
%
%   If P is a two element vector, P(2) is used as a cutoff for signal and
%   noise subspace separation.  All eigenvalues greater than P(2) times
%   the smallest eigenvalue are designated as signal eigenvalues.  In 
%   this case, the signal subspace dimension is at most P(1).
%
%   F = ROOTEIG(...,Fs) uses the sampling frequency Fs in the computation
%   and returns the vector of frequencies, F, in Hz.
%
%   [W,POW] = ROOTEIG(...) returns in addition a vector POW containing the
%   estimates of the powers of the sinusoids in X.
%
%   EXAMPLES:
%      n=0:99;   
%      s=exp(1i*pi/2*n)+2*exp(1i*pi/4*n)+exp(1i*pi/3*n)+randn(1,100);  
%      X=corrmtx(s,12,'mod'); % Estimate the correlation matrix using
%                             % the modified covariance method.
%      [W,P] = rooteig(X,3);     
%   
%   See also ROOTMUSIC, PMUSIC, PEIG, PMTM, PBURG, PWELCH, CORRMTX. 

%   Reference: Stoica, P. and R. Moses, INTRODUCTION TO SPECTRAL ANALYSIS,
%              Prentice-Hall, 1997.

%   Author(s): R. Losada
%   Copyright 1988-2012 The MathWorks, Inc.

narginchk(2,4);

try
   [varargout{1:max(1,nargout)}] = rootmusic(varargin{:},'ev');
catch ME
   throw(ME);
end

% [EOF] rooteig.m

