function [z,p,k] = tf2zpk(b,a)
%TF2ZPK  Discrete-time transfer function to zero-pole conversion.
%   [Z,P,K] = TF2ZPK(NUM,DEN)  finds the zeros, poles, and gain:
%
%                 (z-Z(1))(z-Z(2))...(z-Z(n))
%       H(z) =  K ---------------------------
%                 (z-P(1))(z-P(2))...(z-P(n))
%
%   from a single-input, single-output transfer function in polynomial form:
%
%               NUM(z)
%       H(z) = -------- 
%               DEN(z)
%
%   EXAMPLE:
%     [b,a] = butter(3,.4);
%     [z,p,k] = tf2zpk(b,a)
%
%   See also TF2ZP, ZPLANE.

%   Author(s): T. Bryan
%   Copyright 1988-2008 The MathWorks, Inc.

if nargin < 2 | isempty(a), a = 1; end %#ok

% Must be same length for correct computation of poles/zeros
[b,a] = eqtflength(b,a); 
[z,p,k] = tf2zp(b,a);
  