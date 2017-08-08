function OptimizeScaleValues = setoptimizesv(this, OptimizeScaleValues)
%SETOPTIMIZESV   PreSet function for the 'OptimizeScaleValues' property.

%   Copyright 2008 The MathWorks, Inc.

this.privOptimizeScaleValues = OptimizeScaleValues;

% Quantize the coefficients
quantizecoeffs(this);

% [EOF]
