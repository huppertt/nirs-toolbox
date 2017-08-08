%COEFFS   Returns the filter coefficients in a structure.
%   S = COEFFS(Hd) Returns the filter coefficients of object Hd in the
%   structure S.  The structure will have fields matching the property
%   names in the object Hd.
%
%   If Hd is a multistage object (cascade or parallel), the returned
%   structure will contain fields for each of the stages of the multistage
%   filter.
%
%   % EXAMPLE:
%   b  = fir1(25,.5);
%   Hd = dfilt.dffir(b);
%   c  = coeffs(Hd)
%   b2 = firpm(20,[0 0.4 0.5 1],[1 1 0 0]);
%   Hc = cascade(Hd, dfilt.dffir(b2));
%   c2 = coeffs(Hc)

%   Copyright 1988-2005 The MathWorks, Inc.

% Help for the DFILT method COEFFS.

% [EOF]
