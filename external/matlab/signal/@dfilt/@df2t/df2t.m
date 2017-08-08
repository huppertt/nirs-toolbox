function Hd = df2t(num,den)
%DF2T Direct-Form II Transposed.
%   Hd = DFILT.DF2T(NUM, DEN) constructs a discrete-time direct-form II
%   transposed filter object Hd, with numerator coefficients NUM and
%   denominator coefficients DEN. The leading coefficient of the
%   denominator DEN(1) cannot be 0. 
%
%   Notice that the DSP System Toolbox, along with the Fixed-Point Designer,
%   enables fixed-point support. For more information, see the 
%   <a href="matlab:web([matlabroot,'\toolbox\dsp\dspdemos\html\gsfixedptdemo.html'])">Getting Started with Fixed-Point Filters</a> demo.
%
%   Also, notice that direct-form implementations of IIR filters can lead
%   to numerical problems. In many cases, it can be advantageous to avoid 
%   forming the transfer function and to use a <a href="matlab:help dfilt.df2tsos">second-order section</a>
%   implementation.
%
%   % EXAMPLE #1: Direct instantiation
%   [b,a] = butter(4,.5);
%   Hd = dfilt.df2t(b,a)
%
%   % EXAMPLE #2: Design a 10th order lowpass filter in section order sections
%   f = fdesign.lowpass('N,F3dB',10,.5);  % Specifications
%   Hd = design(f, 'butter', 'Filterstructure', 'df2tsos')
%
%   See also DFILT/STRUCTURES.
  
%   Copyright 1988-2012 The MathWorks, Inc.

Hd = dfilt.df2t;

Hd.FilterStructure = 'Direct-Form II Transposed';

% Hard code the number of coefficients.
Hd.ncoeffs = [1 1];

if nargin>=1
  Hd.Numerator = num;
end
if nargin>=2
  Hd.Denominator = den;
end
