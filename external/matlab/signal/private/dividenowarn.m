function y = dividenowarn(num,den)
% DIVIDENOWARN Divides two polynomials while suppressing warnings.
% DIVIDENOWARN(NUM,DEN) array divides two polynomials but suppresses warnings 
% to avoid "Divide by zero" warnings.

%   Copyright 1988-2002 The MathWorks, Inc.

s = warning; % Cache warning state
warning off  % Avoid "Divide by zero" warnings
y = (num./den);
warning(s);  % Reset warning state

% [EOF] dividenowarn.m