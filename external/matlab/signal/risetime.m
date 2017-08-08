function [r,lwrCross,uprCross,lwrRef,uprRef] = risetime(x, varargin)
%RISETIME Rise time of positive-going bilevel waveform transitions
%   R = RISETIME(X) returns a vector containing the duration of time each
%   transition of the input signal takes to cross from the 10% to 90%
%   reference levels. The signal's sample values are specified by a vector,
%   X, whose sample instants correspond to the index of each sample value.
%   To determine the transitions, the function estimates the state levels
%   of the input waveform by a histogram method and identifies all regions
%   which cross from a 2% tolerance of the lower state level to a 2%
%   tolerance of the upper state level.
% 
%   R = RISETIME(X,Fs) specifies the sample rate, Fs, as a positive
%   scalar, where the first sample instant corresponds to a time of zero.
%
%   R = RISETIME(X,T) specifies the sample instants, T, as a vector with
%   the same number of elements as X.
%
%   R = RISETIME(...,'Tolerance',TOL) specifies the tolerance that the
%   initial and final levels of each positive-going (rising) transition
%   must be within the lower and upper state levels, respectively. The
%   amount, TOL, is a scalar expressed as the percentage of the difference
%   between the upper and lower state levels. The default value is 2.0
%   (percent).
%
%   R = RISETIME(...,'PercentReferenceLevels',PRL) specifies the reference
%   levels as a percentage of the waveform amplitude, where 0% corresponds
%   to the lower state level, and 100% corresponds to the upper state
%   level. PRL is a two-element real row vector whose elements correspond
%   to the lower and upper percent reference levels. The default value for
%   PRL is [10 90].
%
%   R = RISETIME(...,'StateLevels',SL) specifies the levels to use for the
%   lower and upper state levels, SL, as a two-element real row vector
%   whose first and second elements correspond to the lower and upper state
%   levels of the input waveform.
%
%   [R,LT,UT] = RISETIME(...) returns vectors, LT and UT, whose
%   elements correspond to the time instants where X crosses the lower
%   and upper percent reference levels.
%
%   [R,LT,UT,LL,UL] = RISETIME(...) returns the levels, LL and UL, that
%   correspond to the lower and upper percent reference levels.
%
%   RISETIME(...) plots the signal and darkens the regions of each
%   transition where rise time is computed. It marks the lower and upper
%   crossings, and their associated reference levels. The state levels and
%   their associated lower and upper boundaries (adjustable by the TOL
%   parameter) are also plotted.
%
%   % Example 1:
%   %   Plot the rise time information of a 2.3V step waveform 
%   load('transitionex.mat', 'x', 't');
%   risetime(x, t)
%
%   % Example 2:
%   %   Compute the rise time information of a 2.3V step waveform 
%   %   sampled at 4 MHz
%   load('transitionex.mat', 'x', 't');
%   y = risetime(x, 4e6)
%
% See also:  STATELEVELS, FALLTIME, SLEWRATE, OVERSHOOT, MIDCROSS.

%   Copyright 2011-2013 The MathWorks, Inc.

% plot if no output specified
plotFlag = nargout==0;

% obtain durations and crossings for all positive-going transitions
[r,lwrCross,uprCross,lwrRef,uprRef] = transdurs(x,1,plotFlag,varargin{:});


