function Hsos = thissos(Hd, c)
%THISSOS Second-order-section version of this class.
%   Hsos = THISSOS(Hd,C) returns the second-order-section version of the same
%   class as Hd with coefficients in cell array C.  This function is intended
%   as a helper-function to SOS.
%
%   See also DFILT, SOS.

%   Author: Thomas A. Bryan
%   Copyright 1988-2002 The MathWorks, Inc.

Hsos = dfilt.df2sos(c{:});  