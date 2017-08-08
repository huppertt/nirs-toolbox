function [P,DENS] = getNumericSpecs(h)
%GET_SPECS  Get and evaluate design specs from object.

%   Author(s): R. Losada
%   Copyright 1988-2002 The MathWorks, Inc.


Pnorm = get(h,'Pnorm');
initP = get(h,'initPnorm');
P = [initP, Pnorm];
DENS = get(h,'DensityFactor');

