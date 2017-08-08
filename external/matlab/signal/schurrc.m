function [K,E] = schurrc(R)
%SCHURRC Schur algorithm.
%    K = SCHURRC(R) computes the reflection coefficients from autocorrelation
%    vector R. If R is a matrix, SCHURRC finds coefficients for each column of 
%    R, and returns them in the columns of K.
%
%    [K,E] = SCHURRC(R) returns the prediction error variance E. If R is a matrix,
%    SCHURRC finds the error for each column of R, and returns them in the rows of E.
%
%   % Example:
%   %   Create an autocorrelation sequence from the MATLAB speech signal 
%   %   contained in mtlb.mat, and use the Schur algorithm to compute the 
%   %   reflection coefficients of a lattice prediction filter for this 
%   %   autocorrelation sequence.
%
%   load mtlb;                          % Load data
%   r = xcorr(mtlb(1:5),'unbiased');    % Cross correlation
%   k = schurrc(r(5:end))               % Compute reflection coefficients 
%
%   See also LEVINSON. 

%    Copyright 1988-2004 The MathWorks, Inc.
%
%   Reference(s):
% 	  [1] J. Proakis and D. Manolakis, "Digital Signal Processing: Principles,
%         Algorithms, and Applications", pp. 868-873.

error(nargchk(1,1,nargin,'struct'));

% Force column for row vector input:
if (min(size(R)) == 1), R = R(:); end 

[nrows,ncols] = size(R);

for c = 1:ncols,
	A = R(:,c).';	                 % Force row vector
	G = [0 A(2:nrows); A(1:nrows)]; % Initialize generator matrix:

    for m = 2:nrows,
        G(2,:) = [0 G(2,1:nrows-1)];	     % Shift 2nd row of matrix to the right by one
	   K(m-1,c) = -G(1,m)/G(2,m); %#ok       % Compute reflection coefficient
	   G = [1 K(m-1,c); conj(K(m-1,c)) 1]*G; % Update generator matrix
    end
   
	% Return prediction error variance
	E(c,1) = G(2,end);
end


% [EOF] schurrc.m


