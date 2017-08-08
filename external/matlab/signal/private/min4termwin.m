function w = min4termwin(A,x)
%MIN4TERMWIN Generate a minimum 4-term Blackman-Harris window.
%   MIN4TERMWIN(A,X) returns a minimum 4-term Blackman-Harris window using
%   the coefficients specified in A. X is the sampling grid.

%   Author(s): P. Costa 
%   Copyright 1988-2010 The MathWorks, Inc.

%   Reference:
%     [1] fredric j. harris [sic], On the Use of Windows for Harmonic 
%         Analysis with the Discrete Fourier Transform, Proceedings of 
%         the IEEE, Vol. 66, No. 1, January 1978

%w = A(1) - A(2)*cos(x) + A(3)*cos(2.0*x) - A(4)*cos(3.0*x);
B = [A(1); -A(2); A(3); -A(4)];
w = cos(x* [0 1 2 3]) * B;