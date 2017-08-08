function c = coefficientnames(Hd)
%COEFFICIENTNAMES  Coefficient names.
%   COEFFICIENTNAMES(Hd) returns a cell array of the names of the
%   coefficients for this filter structure.
%
%   EXAMPLE:
%           Hd = dfilt.latticemamin;
%           c = coefficientnames(Hd)
%
%   See also DFILT.   
  
%   Author: Thomas A. Bryan
%   Copyright 1988-2008 The MathWorks, Inc.
  
% The singleton filters have extra no-op input parameters so you
% don't have to distinguish the calling syntax between singleton and
% multisection filters for this function.

c = {'Lattice'};
