function c = xcorr2(a,b)
%XCORR2 Two-dimensional cross-correlation.
%   XCORR2(A,B) computes the crosscorrelation of matrices A and B.
%   XCORR2(A) is the autocorrelation function.
%
%   % Example:
%   %   Find the cross-correlation of two matrices a and b:
%   %   a = [2 1 5; 3 1 3; 5 2 2];    b = [1 4 3; 2 5 6];
%   
%   a = [2 1 5; 3 1 3; 5 2 2];    
%   b = [1 4 3; 2 5 6];
%   xcorr2(a,b)
%
%   See also CONV2, XCORR and FILTER2.

%   Author(s): M. Ullman, 2-6-86
%   	   J.N. Little, 6-13-88, revised
%   Copyright 1988-2002 The MathWorks, Inc.

if nargin == 1
	b = a;
end

c = conv2(a, rot90(conj(b),2));

