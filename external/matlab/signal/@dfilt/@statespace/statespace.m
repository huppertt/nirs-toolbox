function Hd = statespace(A,B,C,D)
%STATESPACE Discrete-time, state-space filter.
%   Hd = DFILT.STATESPACE(A, B, C, D) returns a discrete-time state-space
%   filter, Hd, with rectangular arrays A, B, C and D. A, B, C, and D are
%   from the matrix or state-space form of a filter's difference equations:  
%
%   x(n+1) = A*x(n) + B*u(n)
%   y(n)   = C*x(n) + D*u(n)
%
%   where x(n) is the vector states at time n, 
%         u(n) is the input at time n, 
%         y    is the output at time n, 
%         A    is the state-transition matrix, 
%         B    is the input-to-state transmission matrix,
%         C    is the state-to-output transmission matrix, and 
%         D    is the input-to-output transmission matrix. 
%
%   If A, B, C or D are not specified, they default to [], [], [] and 1.
%
%   % EXAMPLE
%   [A,B,C,D] = butter(2,.5);
%   Hd = dfilt.statespace(A,B,C,D)
%
%   See also DFILT/STRUCTURES, TF2SS, ZP2SS   
  
%   Author: Thomas A. Bryan
%   Copyright 1988-2005 The MathWorks, Inc.

Hd = dfilt.statespace;

Hd.FilterStructure = 'State-Space';

% To allow empty inputs to set the defaults
if nargin<1, A = []; end
if nargin<2, B = []; end
if nargin<3, C = []; end
if nargin<4, D = []; end

% Set the defaults if empty input.  This also forces the right sizes of
% empties. 
if isempty(A), A = []; end
if isempty(B), B = zeros(0,1); end
if isempty(C), C = zeros(1,0); end
if isempty(D), D = 1; end

% Validate the consistency of the input before setting the object;
error(abcdchk(A,B,C,D));

Hd.A = A;
Hd.B = B;
Hd.C = C;
Hd.D = D;
