function varargout = pctreflev(stateLevs, varargin)
%PCTREFLEV compute percent reference levels from state levels
%   LEVEL = PCTREFLEV(STATELEVS, PCT) returns an (absolute) reference level
%   that is related to the lower and upper state levels specified by the
%   first and second elements of STATELEVS and a percentage, PCT, where a
%   value of 0 corresponds to the lower state level and a value of 100
%   corresponds to the upper state level:
%
%      REFLEV = STATELEVS(1) + (STATELEVS(2) - STATELEVS(1)) * PCT/100 
% 
%   [REFLEV1, REFLEV2, ...] = PCTREFLEV(STATELEVS, PCT1, PCT2, ...)
%   returns the reference levels that correspond to each percentage
%   specified in the input
%
%   This function is for internal use only. It may be removed in the future.

%   Copyright 2011 The MathWorks, Inc.

%   Reference
%    IEEE Std 181-2003 IEEE Standard on Transitions, Pulses, and Related
%                      Waveforms (Section 5.3.2)

amp = diff(stateLevs)/100;
low = stateLevs(1);
varargout = cell(size(varargin));
for i=1:(nargin-1)
  varargout{i} = low + varargin{i} * amp;
end
