function [b,a] = sos2tf(sos,g)
%SOS2TF 2nd-order sections to transfer function model conversion.
%   [B,A] = SOS2TF(SOS,G) returns the numerator and denominator
%   coefficients B and A of the discrete-time linear system given
%   by the gain G and the matrix SOS in second-order sections form.
%
%   SOS is an L by 6 matrix which contains the coefficients of each
%   second-order section in its rows:
%       SOS = [ b01 b11 b21  1 a11 a21
%               b02 b12 b22  1 a12 a22
%               ...
%               b0L b1L b2L  1 a1L a2L ]
%   The system transfer function is the product of the second-order
%   transfer functions of the sections times the gain G. If G is not
%   specified, it defaults to 1. Each row of the SOS matrix describes
%   a 2nd order transfer function as
%       b0k +  b1k z^-1 +  b2k  z^-2
%       ----------------------------
%       1 +  a1k z^-1 +  a2k  z^-2
%   where k is the row index.
%
%   % Example:
%   %   Compute the transfer function representation of a simple second-
%   %   -order section system.
%
%   sos = [1  1  1  1  0 -1; -2  3  1  1 10  1];
%   [b,a] = sos2tf(sos)
%
%   See also ZP2SOS, SOS2ZP, SOS2SS, SS2SOS

%   Copyright 1988-2013 The MathWorks, Inc.

narginchk(1,2)

if nargin == 1,
    g = 1;
end
% Cast to enforce Precision Rules
if any([signal.internal.sigcheckfloattype(sos,'single','sos2tf','SOS')...
    signal.internal.sigcheckfloattype(g,'single','sos2tf','G')])
  sos = single(sos);
  g = single(g);
end

[L,n] = size(sos);
if n ~= 6,
    error(message('signal:sos2tf:InvalidDimensions'));
end

if L == 0,
    b = [];
    a = [];
    return
end

b = 1;
a = 1;
for m=1:L,
    b1 = sos(m,1:3);
    a1 = sos(m,4:6);
    b = conv(b,b1);
    a = conv(a,a1);
end
b = b.*prod(g);
if length(b) > 3,
    if b(end) == 0,
        b(end) = []; % Remove trailing zeros if any for order > 2
    end
end
if length(a) > 3,
    if a(end) == 0,
        a(end) = []; % Remove trailing zeros if any for order > 2
    end
end

