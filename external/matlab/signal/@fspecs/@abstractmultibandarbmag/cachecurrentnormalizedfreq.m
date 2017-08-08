function cachecurrentnormalizedfreq(this)
%CACHECURRENTNORMALIZEDFREQ

%   Copyright 2011 The MathWorks, Inc.

% Cache the current normalized frequency state before the frequency
% properties are normalized and the privNormalizedFreq property is set to
% true. This will allow us to pass the sampling frequency to the filter
% designer objects that need it to normalize the values specified for the
% BiForcedFrequencyPoints design options.
this.privActualNormalizedFrequencyState = this.privNormalizedFreq;

% [EOF]