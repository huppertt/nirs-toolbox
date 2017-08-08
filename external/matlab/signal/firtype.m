function filtertype = firtype(b)
%FIRTYPE  Type of linear phase FIR filter
%   T = FIRTYPE(B) determines the type, T, (1 through 4) of an FIR filter
%   with coefficients specified in the double or single precision vector B.
%
%   T = FIRTYPE(D) determines the type, T of an FIR digital filter D. You 
%   design a digital filter, D, by calling the <a href="matlab:help designfilt">designfilt</a> function.
%
%   In all cases, the filters must be real and have linear phase.
%
%   Types 1 through 4 are defined as follows:
%
%      Type 1: Even order symmetric coefficients.
%      Type 2: Odd order symmetric coefficients.
%      Type 3: Even order antisymmetric coefficients.
%      Type 4: Odd order antisymmetric coefficients.
%
%   % Example 1:
%   %   Design a 5th-order lowpass FIR filter with normalized cut-off 
%   %   frequency of 0.5 and determine the type of the filter.
%
%   b1 = fir1(5,0.5);    
%   n = firtype(b1)      % odd-order symmetric coefficients
%
%   % Example 2:
%   %   Design a 6th-order lowpass FIR filter with normalized cut-off 
%   %   frequency of 0.5 and determine the type of the filter.
%
%   b2 = fir1(6,0.5);    
%   n = firtype(b2)      % even-order symmetric coefficients
%
%   % Example 3:
%   %   Use the designfilt function to design a 10-th order FIR Hilbert
%   %   filter with a transition width of 0.1*pi radians/s. Verify that the
%   %   FIR type is 3. 
%
%   D = designfilt('hilbertfir', 'FilterOrder', 10, 'TransitionWidth', 0.1);
%   n = firtype(D)     
%
%   See also FVTOOL, FILTORD, ISALLPASS, ISLINPHASE, ISMAXPHASE, ISMINPHASE

%   Copyright 2012-2013 The MathWorks, Inc.

narginchk(1,1);

if all(size(b)>[1 1])
  error(message('signal:signalanalysisbase:inputnotsupportedFIR'));
end

if ~islinphase(b),
  error(message('signal:signalanalysisbase:expectLinearPhase'));
end

if ~isreal(b),
  error(message('signal:signalanalysisbase:mustBeReal'));
end

h = impz(b);

if length(h) == 1,
  filtertype = 1;
else
  if strcmpi('symmetric',signalpolyutils('symmetrytest',h,1)),
    % Symmetric coefficients
    issymflag = 1;
  else
    % We already know it is linear phase, must be antisymmetric
    % coefficients
    issymflag = 0;
  end  
  filtertype = signalpolyutils('determinetype',h,issymflag,1);
end


