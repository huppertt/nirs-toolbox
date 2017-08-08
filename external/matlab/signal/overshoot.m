function [os,osLev,osInst] = overshoot(x, varargin)
%OVERSHOOT Overshoot metrics of bilevel waveform transitions
%   OS = OVERSHOOT(X) returns the signal's greatest deviation larger than
%   the final state level of each transition. The overshoots are expressed
%   as a percentage of the difference between the state levels as a vector,
%   OS, whose length corresponds to the number of transitions detected in
%   the input signal. The signal's sample values are specified by a vector,
%   X, whose sample instants correspond to the index of each sample value.
%   To determine the transitions, the function estimates the state levels
%   of the input waveform by a histogram method and identifies all regions
%   which cross between a 2% tolerance of each state level.
%
%   OS = OVERSHOOT(X,Fs) specifies the sample rate, Fs, as a positive
%   scalar, where the first sample instant corresponds to a time of zero.
%
%   OS = OVERSHOOT(X,T) specifies the sample instants, T, as a vector with
%   the same number of elements as X.
%
%   OS = OVERSHOOT(...,'Region',REGION) specifies the region, REGION, over
%   which to perform the overshoot computation as one of 'Preshoot' |
%   'Postshoot', where the default value is 'Postshoot'. If 'Preshoot' is
%   specified, the end of the pre-transition aberration region of each is
%   defined as the last instant where the signal exits the first state. If
%   'Postshoot' is specified, the start of the post-transition aberration
%   region is defined as the instant when the signal enters the second
%   state.
%
%   OS = OVERSHOOT(...,'SeekFactor',FACTOR) specifies the duration of the
%   region over which to compute the overshoot for each transition as a
%   multiple of the corresponding transition duration. If the edge of the
%   waveform is reached or a complete intervening transition is detected
%   before this duration elapses, the duration is truncated to the edge of
%   the waveform or the start of the intervening transition. The default
%   value of this parameter is 3.
% 
%   OS = OVERSHOOT(...,'Tolerance',TOL) specifies the tolerance that the
%   initial and final levels of each transition must be within their
%   respective state levels. The amount, TOL, is a scalar expressed
%   as the percentage of the difference between the upper and lower state
%   levels. The default value is 2.0 (percent).
%
%   OS = OVERSHOOT(...,'PercentReferenceLevels',PRL) specifies the
%   reference levels as a percentage of the waveform amplitude, where 0%
%   corresponds to the lower state level, and 100% corresponds to the upper
%   state level. PRL is a two-element real row vector whose elements
%   correspond to the lower and upper percent reference levels. The default
%   value of this parameter is [10 90].
% 
%   OS = OVERSHOOT(...,'StateLevels',SL) specifies the levels to use for
%   the lower and upper state levels, SL, as a two-element real row vector
%   whose first and second elements correspond to the lower and upper state
%   levels of the input waveform.
% 
%   [OS,OSLEV,OSINST] = OVERSHOOT(...) returns vectors, OSLEV and OSINST,
%   whose elements correspond to the level and sample instant of the
%   overshoot of each transition.
%
%   OVERSHOOT(...) plots the signal and marks the location of the overshoot
%   of each transition as well as the lower and upper crossings and their
%   associated reference levels. The state levels and their associated
%   lower and upper boundaries (adjustable by the TOL parameter) are also
%   plotted.
%
%   % Example 1:
%   %   Plot overshoot information of a 2.3V step waveform 
%   load('transitionex.mat', 'x', 't');
%   overshoot(x, t)
%
%   % Example 2:
%   %   Compute overshoot information of a 2.3V step waveform 
%   %   sampled at 4 MHz
%   load('transitionex.mat', 'x', 't');
%   y = overshoot(x, 4e6)
%
% See also:  STATELEVELS, UNDERSHOOT, RISETIME, SETTLINGTIME, MIDCROSS.

%   Copyright 2011-2013 The MathWorks, Inc.

% plot if no output specified
plotFlag = nargout==0;

% obtain overshoot information
[os,osLev,osInst] = transshoot(x,1,plotFlag,varargin{:});