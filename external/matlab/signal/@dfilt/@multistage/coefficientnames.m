function c = coefficientnames(Hd)
%COEFFICIENTNAMES  Coefficient names.
%   COEFFICIENTNAMES(Hd) returns a vector of cell arrays of the names of the
%   coefficients each section of the discrete-time filter Hd.
%   Example:
%     Hd = cascade(dfilt.df2t,dfilt.latticear);
%     c  = coefficientnames(Hd)
%
%   See also DFILT.   
  
%   Copyright 1988-2014 The MathWorks, Inc.
  
c = cell(1,length(Hd.Stage));
for k=1:length(Hd.Stage)
  c{k} = coefficientnames(Hd.Stage(k));
end

