function Href = reffilter(this)
%REFFILTER   Reference double-precision floating-point filter.
%   Href = REFFILTER(Hd) returns a new filter Href that has the same
%   structure as Hd, but uses the reference coefficients and has its
%   arithmetic property set to 'double'.
%
%   REFFILTER(Hd) differs from DOUBLE(Hd) in that the returned filter has
%   the reference coefficients of Hd instead of the actual (quantized)
%   coefficients.
%
%   REFFILTER(Hd) is useful to quickly have a reference version of a
%   filter one is working with for comparison purposes.
%
%   % Example: Compare several fixed-point filters with double-precision
%   % floating-point filter.
%   f = fdesign.lowpass('N,Fc,Ap,Ast',87,.5,1e-3,1e-6,'linear');
%   H = design(f,'FilterStructure','dffir');
%   H1 = copy(H); H2 = copy(H); % Generate two copies of the same filter
%   H.Arithmetic = 'fixed';  % Set H to filter using fixed-point arithmetic
%   H1.Arithmetic = 'fixed'; % Same for H1
%   H2.Arithmetic = 'fixed'; % Same for H2
%   H.CoeffWordLength  = 16;  % Use 16 bits to represent the coefficients
%   H1.CoeffWordLength = 12;  % Use 12 bits to represent the coefficients
%   H2.CoeffWordLength =  8;  % Use  8 bits to represent the coefficients
%   Href = reffilter(H);
%   Hfvt = fvtool(Href,H,H1,H2); 
%   set(Hfvt,'ShowReference','off'); % Reference already displayed (once)
%   legend(Hfvt,'Reference filter','16-bits','12-bits','8-bits');
%
%   See also dfilt/double.

%   Author(s): R. Losada
%   Copyright 2003-2005 The MathWorks, Inc.



% [EOF]
