function [anxt,enxt]=levup(acur,knxt,ecur)
%LEVUP  One step forward Levinson recursion
%   Anxt=LEVUP(Acur,Knxt) returns the P+1'th order prediction polynomial,
%   Anxt based on the P'th order prediction polynomial, Acur, and the 
%   P+1'th order reflection coefficient, Knxt.
%
%   [Anxt,Enxt]=LEVUP(Acur,Knxt,Ecur) returns the P+1'th order prediction
%   prediction error, Enxt based on the P'th order prediction error, Ecur.  
%
%   See also LEVDOWN, LEVINSON, RLEVINSON.

%   References:  P. Stoica R. Moses, Introduction to Spectral Analysis
%                Prentice Hall, N.J., 1997, Chapter 3.
%
%   Author(s): A. Ramasubramanian,
%   Copyright 1988-2002 The MathWorks, Inc.

% Some preliminaries first


if (nargout==2 & nargin<3) 
   error(message('signal:levup:Nargchk'));
end

% To convert to a column vector if not already so
acur = acur(:);                  
acur = acur(2:end);  % Drop the leading 1, it is not needed in the stepup

% Matrix formulation from Stoica is used to avoid looping
anxt = [acur;0]+knxt*conj([acur(end:-1:1);1]);

if nargin==3
  enxt = (1-knxt'.*knxt)*ecur;
end

anxt = [1;anxt];     % Append the one to make it a true polynomial
anxt = anxt.';       % Return the polynomial as a row vector

% [EOF] levup.m
