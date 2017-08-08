function this = allpass(c)
%ALLPASS Minimum-multiplier allpass filter.
%   Hd = DFILT.ALLPASS(C) constructs a minimum-multiplier allpass structure
%   given the allpass coefficients in vector C.
%
%   C must have between one and four coefficients. 
%
%   The allpass transfer function of Hd given the coefficients in C is:
%                           -1          -n 
%             C(n) + C(n-1)z + .... +  z
%     H(z) = -------------------------------
%                      -1             -n
%             1 + C(1)z + .... + C(n)z
%
%   Notice that the leading '1' coefficient in the denominator is not
%   part of C.
%
%   It is possible to construct a cascade of these filters using
%   DFILT.CASCADEALLPASS. See the help for that filter structure for more
%   information.
%
%   Example: Construct a second-order minimum-multiplier allpass filter
%   C = [1.5,0.7];
%   Hd = dfilt.allpass(C);
%   info(Hd)
%   realizemdl(Hd) % Requires Simulink; build model for filter
%
%   See also DFILT/STRUCTURES

%   Author(s): R. Losada
%   Copyright 2005-2006 The MathWorks, Inc.

this = dfilt.allpass;

this.FilterStructure = 'Minimum-Multiplier Allpass';

if nargin > 0,
    if ~isreal(c),
        error(message('signal:dfilt:allpass:allpass:complexCoeffs'));
    end
    this.AllpassCoefficients = c;
end

% [EOF]
