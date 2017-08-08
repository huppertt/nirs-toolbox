function Hd = latticemamin(lattice)
%LATTICEMAMIN Lattice Moving-Average for Minimum Phase.
%   Hd = DFILT.LATTICEMAMIN(K) constructs a Lattice moving-average (MA) for
%   minimum phase discrete-time filter object Hd with lattice coefficients K.
%   If K is not specified, it defaults to []. In this case, the filter
%   passes the input through to the output unchanged.
%
%   Notice that if the K coefficients define a minimum phase filter, the
%   resulting filter in this structure is minimum phase. If your
%   coefficients do not define a minimum phase filter, placing them in this
%   structure does not produce a minimum phase filter. 
%
%   Also, the DSP System Toolbox, along with the Fixed-Point Designer,
%   enables fixed-point support. For more information, see the 
%   <a href="matlab:web([matlabroot,'\toolbox\dsp\dspdemos\html\gsfixedptdemo.html'])">Getting Started with Fixed-Point Filters</a> demo.
%
%   % EXAMPLE
%   k = [.66 .7 0.44];
%   Hd = dfilt.latticemamin(k)
%   ismin = isminphase(Hd)
%   realizemdl(Hd); % Requires Simulink
%
%   See also DFILT/STRUCTURES, TF2LATC
  
%   Copyright 1988-2012 The MathWorks, Inc.
Hd = dfilt.latticemamin;

Hd.FilterStructure = 'Lattice Moving-Average (MA) For Minimum Phase';

if nargin>=1
  Hd.Lattice = lattice;
end
