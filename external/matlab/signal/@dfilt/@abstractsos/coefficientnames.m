function c = coefficientnames(Hd)
%COEFFICIENTNAMES  Coefficient names.
%   COEFFICIENTNAMES(Hd) returns a cell array of the names of the
%   coefficients for this filter structure.
%
%   Example:
%     Hd = dfilt.df1sos;
%     c = coefficientnames(Hd)
%
%   See also DFILT.   
  
%   Author: R. Losada
%   Copyright 1988-2003 The MathWorks, Inc.
  
% The singleton filters have extra no-op input parameters so you
% don't have to distinguish the calling syntax between singleton and
% multisection filters for this function.

c = {'SOSMatrix', 'ScaleValues'};