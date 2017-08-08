%FIRTYPE  Determine the type of a linear phase FIR filter.
%   T = FIRTYPE(Hd) determines the type (1 through 4) of a Discrete-Time
%   FIR filter object Hd.  The filter must be real and have linear phase.
%
%   Type 1 through 4 are defined as follows:
%      
%     - Type 1: Even order symmetric coefficients.
%     - Type 2: Odd order symmetric coefficients.
%     - Type 3: Even order antisymmetric coefficients.
%     - Type 4: Odd order antisymmetric coefficients.
%
%   If Hd is a cascade or parallel filter, all sections must be real FIR       
%   filters with linear phase.  For these cases, T is a cell array
%   containing the type of each section. 
%
%   See also DFILT.

% Copyright 1988-2004 The MathWorks, Inc.

% Help for the filter's FIRTYPE method.

% [EOF]
