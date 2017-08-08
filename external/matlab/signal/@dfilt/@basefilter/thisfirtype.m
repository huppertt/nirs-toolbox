function filtertype = thisfirtype(Hd)
%FIRTYPE  Determine the type (1-4) of a linear phase FIR filter.
%   T = FIRTYPE(Hd) determines the type (1 through 4) of an FIR filter
%   object Hd.  The filter must be real and have linear phase.
%
%   Type 1 through 4 are defines as follows:
%      
%     - Type 1: Even order symmetric coefficients.
%     - Type 2: Odd order symmetric coefficients.
%     - Type 3: Even order antisymmetric coefficients.
%     - Type 4: Odd order antisymmetric coefficients.
%
%   If Hd has multiple sections, all sections must be real FIR filters with
%   linear phase.  For this case, T is a cell array containing the type
%   of each section.

%   Author(s): R. Losada, T. Bryan
%   Copyright 1988-2004 The MathWorks, Inc.

if ~isfir(Hd),
    error(message('signal:dfilt:basefilter:thisfirtype:expectFIR'));
end

if ~islinphase(Hd),
    error(message('signal:dfilt:basefilter:thisfirtype:expectLinearPhase'));
end

if ~isreal(Hd),
    error(message('signal:dfilt:basefilter:thisfirtype:MustBeReal'));
end

% The FIR type is based on the impulse response.  For direct-form FIR
% filters, this is just the numerator coefficients.  However, for other FIR
% filters, such as lattice-ma, the impulse response has a complicated
% relationship to the coefficients, and so impz is the best for all cases.
h = impz(Hd);

if length(h) == 1,
  filtertype = 1;
else
    if strcmpi('symmetric',signalpolyutils('symmetrytest',h,1)),
        % Symmetric coefficients
        issymflag = 1;
    else
        % We already know it is linear phase, must be antisymmetric coefficients
        issymflag = 0;
    end
    
    filtertype = signalpolyutils('determinetype',h,issymflag,1);
end

% [EOF]
