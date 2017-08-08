function Hd = df1(num,den)
%DF1 Direct-Form I.
%   Hd = DFILT.DF1(NUM, DEN) constructs a discrete-time direct-form I
%   filter object Hd, with numerator coefficients NUM and denominator
%   coefficients DEN. The leading coefficient of the denominator DEN(1)
%   cannot be 0.
%
%   Notice that the DSP System Toolbox, along with the Fixed-Point Designer,
%   enables fixed-point support. For more information, see the 
%   <a href="matlab:web([matlabroot,'\toolbox\dsp\dspdemos\html\gsfixedptdemo.html'])">Getting Started with Fixed-Point Filters</a> demo.
%
%   Also, notice that direct-form implementations of IIR filters can lead
%   to numerical problems. In many cases, it can be advantageous to avoid 
%   forming the transfer function and to use a <a href="matlab:help dfilt.df1sos">second-order section</a>
%   implementation.
%
%   % EXAMPLE #1: Direct instantiation
%   [b,a] = butter(4,.5);
%   Hd = dfilt.df1(b,a)
%
%   % EXAMPLE #2: Design a 10th order lowpass filter in section order sections
%   f = fdesign.lowpass('N,F3dB',10,.5);  % Specifications
%   Hd = design(f, 'butter', 'Filterstructure', 'df1sos')
%
%   See also DFILT/STRUCTURES.
  
%   Copyright 1988-2012 The MathWorks, Inc.

Hd = dfilt.df1;
Hd.ncoeffs = [1 1];

Hd.FilterStructure = 'Direct-Form I';

% Tap Index is a vector of two elements. The first element corresponds to 
% the WRITE index for the numerator circular buffer and the second element 
% corresponds to the WRITE index for the denominator circular buffer.
Hd.tapIndex = [0 0];

% Hard code the number of coefficients to avoid special cases in the 
% thissetstates and getstates methods.
Hd.ncoeffs = [1 1];

if nargin>=1
  Hd.Numerator = num;
end
if nargin>=2
  Hd.Denominator = den;
end
