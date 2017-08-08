function filtertype = firtype(Hb)
%FIRTYPE  Determine the type (1-4) of a linear phase FIR filter.
%   T = FIRTYPE(Hb) determines the type (1 through 4) of an FIR filter
%   object Hb.  The filter must be real and have linear phase.
%
%   Type 1 through 4 are defines as follows:
%      
%     - Type 1: Even order symmetric coefficients.
%     - Type 2: Odd order symmetric coefficients.
%     - Type 3: Even order antisymmetric coefficients.
%     - Type 4: Odd order antisymmetric coefficients.
%
%   If Hb has multiple sections, all sections must be real FIR filters with
%   linear phase.  For this case, T is a cell array containing the type
%   of each section.

%   Author(s): R. Losada, T. Bryan
%   Copyright 1988-2003 The MathWorks, Inc.

filtertype = base_num(Hb, 'thisfirtype');
