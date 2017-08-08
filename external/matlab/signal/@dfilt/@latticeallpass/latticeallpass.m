function Hd = latticeallpass(lattice)
%LATTICEALLPASS Lattice Allpass.
%   Hd = DFILT.LATTICEALLPASS(LATTICE) constructs a discrete-time lattice 
%   allpass filter object Hd with lattice coefficients K. If K is not
%   specified, it defaults to []. In this case, the filter passes the input
%   through to the output unchanged. 
%
%   Notice that the DSP System Toolbox, along with the Fixed-Point Designer,
%   enables fixed-point support. For more information, see the 
%   <a href="matlab:web([matlabroot,'\toolbox\dsp\dspdemos\html\gsfixedptdemo.html'])">Getting Started with Fixed-Point Filters</a> demo.
%
%   % EXAMPLE
%   k = [.66 .7 .44];
%   Hd = dfilt.latticeallpass(k)
%   realizemdl(Hd); % Requires Simulink
%
%   See also DFILT/STRUCTURES, TF2LATC
  
%   Copyright 1988-2012 The MathWorks, Inc.

Hd = dfilt.latticeallpass;

Hd.FilterStructure = 'Lattice Allpass';

if nargin>=1
  Hd.Lattice = lattice;
end
