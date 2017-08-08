%SOS  Convert to second-order-sections (for IIR filters only).
%   Hsos = SOS(Hd) converts IIR discrete-time filter Hd to second-order section
%   form. Hd can be of class dfilt.df1, dfilt.df1t, dfilt.df2 or
%   dfilt.df2t.
%
%   SOS(Hd,DIR_FLAG) specifies the ordering of the 2nd order sections. If
%   DIR_FLAG is equal to 'UP', the first row will contain the poles closest
%   to the origin, and the last row will contain the poles closest to the
%   unit circle. If DIR_FLAG is equal to 'DOWN', the sections are ordered
%   in the opposite direction. The zeros are always paired with the poles
%   closest to them. DIR_FLAG defaults to 'UP'.
%
%     % Example 1:
%     [b,a] = butter(8,.5);
%     Hd = dfilt.df2(b,a);
%     Hsos = sos(Hd,'up',inf)
%   
%     % Example 2: 
%     % It is preferable to design directly in SOS form rather than to
%     % construct the transfer function and then convert which can
%     % introduce roundoff error
%     f = fdesign.lowpass('N,F3db',8,0.5);
%     Hsos = design(f,'butter','FilterStructure','df1sos');

%   Copyright 1988-2010 The MathWorks, Inc.

% Help for the p-coded SOS method of DFILT.DF1, DFILT.DF1T, DFILT.DF2 and
% DFILT.DF2T classes.
