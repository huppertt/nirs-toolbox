function y=diric(x,N)
%DIRIC	Dirichlet, or periodic sinc function
%   Y = DIRIC(X,N) returns a matrix the same size as X whose elements
%   are the Dirichlet function of the elements of X.  Positive integer
%   N is the number of equally spaced extrema of the function in the 
%   interval 0 to 2*pi. 
%
%   The Dirichlet function is defined as
%       d(x) = sin(N*x/2)./(N*sin(x/2))   for x not a multiple of 2*pi
%              +1 or -1 for x a multiple of 2*pi. (depending on limit)
% 
%   % Example 1:
%   %   Plot the Dirichlet function over the range 0 to 4? for N = 7 and 
%   %   N = 8. 
%
%   x = linspace(0,4*pi,300);           % Generate linearly spaced vectors
%   subplot(211); plot(x,diric(x,7)); 
%   title('Diric, N = 7'); axis tight;     
%   subplot(212); plot(x,diric(x,8));      
%   title('Diric, N = 8'); axis tight; 
%
%   % Example 2:
%   %   Plot and display the difference in shape of Dirichlet and Sinc 
%   %   functions. 
%
%   x_diric = linspace(0,4*pi,300);     % Generate linearly spaced vectors
%   x_sinc  = linspace(-5,5);           % Generate linearly spaced vectors
%   subplot(211); plot(x_diric,diric(x_diric,7));
%   title('Diric, N = 7'); axis tight;     
%   subplot(212); plot(x_sinc,sinc(x_sinc)); 
%   title('Sinc, Range: -5 to 5 '); axis tight;         

%   Author(s): T. Krauss, 1-14-93
%   Copyright 1988-2004 The MathWorks, Inc.

error(nargchk(2,2,nargin,'struct'));
if round(N) ~= N || N < 1 || numel(N) ~= 1,
    error(message('signal:diric:MustBePosInteger'));
end

y=sin(.5*x);
i=find(abs(y)>1e-12);   % set where x is not divisible by 2 pi
j=1:length(x(:));
j(i)=[];                         % complement set
y(i)=sin((N/2)*x(i))./(N*y(i));
y(j)=sign(cos(x(j)*((N+1)/2)));

