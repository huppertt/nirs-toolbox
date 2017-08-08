function [targs, strs] = convertchoices(this)
%CONVERTCHOICES

%   Copyright 2012 The MathWorks, Inc.

% Common filter structures for both dfilt 
strs = {'Direct-Form I',...
    'Direct-Form II',...
    'Direct-Form I Transposed',...
    'Direct-Form II Transposed',...
    'State-Space',...
    'Lattice Autoregressive Moving-Average (ARMA)'};

targs = {'df1','df2','df1t','df2t','statespace','latticearma'};

refHd = reffilter(this);

% FD Tbx required for coupled-allpass conversions
strs = [strs ,...
  {'Coupled-Allpass (CA) Lattice',...
  'Coupled-Allpass (CA) Lattice with Power-Complementary (PC) Output'}];
targs = [targs {'calattice', 'calatticepc'}];

if isallpass(refHd),
  strs = [strs {'Lattice allpass'}];
  targs = [targs {'latticeallpass'}];
end


% [EOF]
