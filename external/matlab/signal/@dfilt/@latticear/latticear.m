function Hd = latticear(lattice)
%LATTICEAR Lattice Autoregressive (AR).
%   Hd = DFILT.LATTICEAR(LATTICE) constructs a discrete-time lattice AR
%   filter object Hd with lattice coefficients K. If K is not
%   specified, it defaults to []. In this case, the filter passes the input
%   through to the output unchanged. 
%
%   Notice that the DSP System Toolbox, along with the Fixed-Point Designer,
%   enables fixed-point support. For more information, see the 
%   <a href="matlab:web([matlabroot,'\toolbox\dsp\dspdemos\html\gsfixedptdemo.html'])">Getting Started with Fixed-Point Filters</a> demo.
%
%   % EXAMPLE
%   k = [.66 .7 .44];
%   Hd = dfilt.latticear(k)
%   realizemdl(Hd); % Requires Simulink
%
%   See also DFILT/STRUCTURES, TF2LATC
  
%   Copyright 1988-2012 The MathWorks, Inc.

Hd = dfilt.latticear;

Hd.FilterStructure = 'Lattice Autoregressive (AR)';

if nargin>=1
  Hd.Lattice = lattice;
end
