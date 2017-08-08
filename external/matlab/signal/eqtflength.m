function [b,a,N,M] = eqtflength(b,a)
%EQTFLENGTH   Equalize the length of a discrete-time transfer function.
%   [B,A] = EQTFLENGTH(NUM,DEN) forces NUM and DEN to be of the same
%   length by appending zeros to either one as necessary.  If both NUM
%   and DEN have common trailing zeros, they are removed from both of
%   them.
%
%   EQTFLENGTH is intended to be used with discrete-time transfer
%   functions expressed in terms of negative powers of Z only.
%
%   [B,A,N,M] = EQTFLENGTH(NUM,DEN) returns the numerator order N and
%   the denominator order M, not including trailing zeros.

%   Author(s): R. Losada
%   Copyright 1988-2004 The MathWorks, Inc.

error(nargchk(2,2,nargin,'struct'));

% First make b and a rows
b = b(:).';
a = a(:).';

% Catch cases when the den is zero or empty.
if isempty(a) || max(abs(a)) == 0,
    % Divide by zero not allowed
    error(message('signal:eqtflength:InvalidRange'));
end

% First make them of equal length
a = [a zeros(1,max(0,length(b)-length(a)))];
b = [b zeros(1,max(0,length(a)-length(b)))];

% Now remove trailing zeros, but only if present in both b and a
i = find(a~=0);
j = find(b~=0);

% Get the orders of the numerator and denominator
M = i(end) - 1;

% If the numerator is all zeros, j will be empty, catch this case
% note that i will never be empty, if a is all zeros, an error is 
% returned above.
if isempty(j),
	N = 0;
else
	N = j(end) - 1;
end

% Get the index of the largest negative order nonzero element
n = max(M+1,N+1);

a = a(1:n);
b = b(1:n);

% [EOF] - eqtflength.m

