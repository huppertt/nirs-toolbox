function h = dfiir(numstates,denstates)
%DFIIR   Direct-form IIR filter states.
%   H = FILTSTATES.DFIIR constructs a default direct-form IIR filter states
%   object.
%
%   H = FILTSTATES.DFIIR(NUMSTATES,DENSTATES) constructs an object and sets
%   its 'Numerator' and 'Denominator' properties to NUMSTATES and DENSTATES
%   respectively.  
%
%   Notice that the DSP System Toolbox, along with the Fixed-Point
%   Designer, enables single precision floating-point and fixed-point
%   support for the Numerator and Denominator states.
%
%   Example #1, construct the default object
%   h = filtstates.dfiir
%
%   Example #2, construct an object with Numerator and Denominator states
%   as vectors of zeros.
%   h = filtstates.dfiir(zeros(4,1),zeros(4,1));
%
%   See also FILTSTATES.DOUBLE, DFILT.

%   Author(s): P. Costa
%   Copyright 1988-2012 The MathWorks, Inc.

error(nargchk(0, 2, nargin,'struct'));

h = filtstates.dfiir;

if nargin>=1
  h.Numerator = numstates;
end
if nargin>=2
  h.Denominator = denstates;
end

% [EOF]
