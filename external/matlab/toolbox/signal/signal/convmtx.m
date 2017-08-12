function t = convmtx(v,n)
%CONVMTX Convolution matrix.
%   CONVMTX(C,N) returns the convolution matrix for vector C. If C is a
%   column vector and X is a column vector of length N, then CONVMTX(C,N)*X
%   is the same as CONV(C,X). If R is a row vector and X is a row vector of
%   length N, then X*CONVMTX(R,N) is the same as CONV(R,X).
%
%   % Example:
%   %   Generate a simple convolution matrix.
%
%   h = [1 2 3 2 1];
%   convmtx(h,7)        % Convolution matrix
%
%   See also CONV.

%   Copyright 1988-2013 The MathWorks, Inc.

% Cast to enforce Precision Rules
n = signal.internal.sigcasttofloat(n,'double','convmtx','N','allownumeric');

% Checks if 'C' is a valid numeric data input
isVsingle = signal.internal.sigcheckfloattype(v,'single','convmtx','C',...
  'allownumeric');

[mv,nv] = size(v);
v = v(:); % make v a column vector

%  t = toeplitz([v; zeros(n-1,1)],zeros(n,1));  put Toeplitz code inline
c = [v; zeros(n-1,1)]; 
r = zeros(n,1);
m = length(c);
x = [r(n:-1:2) ; c(:)]; % build vector of user data
%
cidx = (0:m-1)';
ridx = n:-1:1;
t = cidx(:,ones(n,1)) + ridx(ones(m,1),:); % Toeplitz subscripts
t(:) = x(t); % actual data
if isVsingle
  t = single(t);
end
% end of toeplitz code

if mv < nv
    t = t.';
end

