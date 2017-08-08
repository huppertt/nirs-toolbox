function [acur,ecur,knxt]=levdown(anxt,enxt)
%LEVDOWN  One step backward Levinson recursion
%   Acur=LEVDOWN(Anxt) returns the P'th order prediction polynomial,                
%   Acur based on the P+1'th order prediction polynomial, Anxt.  
%
%   [Acur,Ecur]=LEVDOWN(Anxt,Enxt) returns the P'th order prediction 
%   error, Ecur, based on the P+1'th order prediction error, Enxt. 
%
%   [Acur,Ecur,Knxt]=LEVDOWN(...) returns the P+1'th order reflection
%   coefficient, Knxt.
%
%
%   See also LEVUP, LEVINSON, RLEVINSON.

%   References: P. Stoica R. Moses, Introduction to Spectral Analysis
%               Prentice Hall, N.J., 1997, Chapter 3.
%
%   Author(s): A. Ramasubramanian 
%   Copyright 1988-2002 The MathWorks, Inc.

% Some preliminaries first
if nargout>=2 & nargin<2
   error(message('signal:levdown:Nargchk'));
end

% Convert to a column vector if not already so
anxt = anxt(:);                  
anxt = anxt(2:end);   % Drop the leading 1, it is not needed 
                      % in the step down

% Extract the k+1'th reflection coefficient
knxt = anxt(end);
if knxt == 1.0,
   error(message('signal:levdown:InvalidParam'));
end

% A Matrix formulation from Stoica is used to avoid looping 
acur = (anxt(1:end-1)-knxt*conj(anxt(end-1:-1:1)))/(1-abs(knxt)^2);

if nargin == 2
   ecur = enxt/(1-knxt'.*knxt);
end

acur = [1;acur];      % Append the one to make it a true polynomial
acur = acur.';        % Return the polynomial as a row vector

% [EOF] levdown.m
