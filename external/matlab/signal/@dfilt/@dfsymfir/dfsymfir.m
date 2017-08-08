function Hd = dfsymfir(num)
%DFSYMFIR Direct-Form Symmetric FIR.
%   Hd = DFILT.DFSYMFIR(NUM) constructs a discrete-time, direct-form
%   symmetric FIR filter object Hd, with numerator coefficients NUM. 
%  
%   Note that one usually does not construct DFILT filters explicitly.
%   Instead, one obtains these filters as a result from a design using <a
%   href="matlab:help fdesign">FDESIGN</a>. 
%
%   Also, the DSP System Toolbox, along with the Fixed-Point Designer,
%   enables fixed-point support. For more information, see the 
%   <a href="matlab:web([matlabroot,'\toolbox\dsp\dspdemos\html\gsfixedptdemo.html'])">Getting Started with Fixed-Point Filters</a> demo.
%
%   % EXAMPLE #1: Direct instantiation
%   b = [-0.008 0.06 0.44 0.44 0.06 -0.008];
%   Hd = dfilt.dfsymfir(b)
%   realizemdl(Hd)    % Requires Simulink
%   
%   % EXAMPLE #2: Design an equiripple lowpass filter with default specifications
%   Hd = design(fdesign.lowpass, 'equiripple', 'Filterstructure', 'dfsymfir');
%   fvtool(Hd)        % Analyze filter
%   x = randn(100,1); % Input signal
%   y = filter(Hd,x); % Apply filter to input signal
%
%   See also DFILT/STRUCTURES
  
%   Copyright 1988-2012 The MathWorks, Inc.

Hd = dfilt.dfsymfir;
Hd.ncoeffs = 1;
Hd.HiddenStates = 0;
Hd.TapIndex = 0;

Hd.FilterStructure = 'Direct-Form Symmetric FIR';

if nargin>=1
  Hd.Numerator = num;
end
