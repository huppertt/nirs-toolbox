function [f,lwrCross,uprCross,lwrRef,uprRef] = falltime(x, varargin)
%FALLTIME Fall time of negative-going bilevel waveform transitions
%   F = FALLTIME(X) returns a vector containing the duration of time each
%   transition of the input signal takes to cross from the 90% to 10%
%   reference levels. The signal's sample values are specified by a vector,
%   X, whose sample instants correspond to the index of each sample value.
%   To determine the transitions, the function estimates the state levels
%   of the input waveform by a histogram method and identifies all regions
%   which cross from a 2% tolerance of the upper state level to a 2%
%   tolerance of the lower state level.
% 
%   F = FALLTIME(X,Fs) specifies the sample rate, Fs, as a positive scalar,
%   where the first sample instant corresponds to a time of zero.
%
%   F = FALLTIME(X,T) specifies the sample instants, T, as a vector with the
%   same number of elements as X. 
%
%   F = FALLTIME(...,'Tolerance',TOL) specifies the tolerance that the
%   initial and final levels of each negative-going (falling) transition
%   must be within the upper and lower state levels, respectively. The
%   amount, TOL, is a scalar expressed as the percentage of the difference
%   between the upper and lower state levels. The default value is 2.0
%   (percent).
%
%   F = FALLTIME(...,'PercentReferenceLevels',PRL) specifies the reference
%   levels as a percentage of the waveform amplitude, where 0% corresponds
%   to the lower state level, and 100% corresponds to the upper state
%   level. PRL is a two-element real row vector whose elements correspond
%   to the lower and upper percent reference levels. The default value for
%   PRL is [10 90].
%
%   F = FALLTIME(...,'StateLevels',SL) specifies the levels to use for the
%   lower and upper state levels, SL, as a two-element real row vector
%   whose first and second elements correspond to the lower and upper state
%   levels of the input waveform.
%
%   [F,LT,UT] = FALLTIME(...) returns vectors, LT and UT, whose elements
%   correspond to the time instants where X crosses the lower and upper
%   percent reference levels.
%
%   [F,LT,UT,LL,UL] = FALLTIME(...) returns the levels, LL and UL, that
%   correspond to the lower and upper percent reference levels.
%
%   FALLTIME(...) plots the signal and darkens the regions of each
%   transition where fall time is computed. It marks the lower and upper
%   crossings, and their associated reference levels. The state levels and
%   their associated lower and upper boundaries (adjustable by the TOL
%   parameter) are also plotted.
%
%   % Example 1:
%   %   Plot the fall time information of a 2.3V step waveform 
%   load('negtransitionex.mat', 'x', 't');
%   falltime(x, t)
%
%   % Example 2:
%   %   Compute the fall time information of a 2.3V step waveform 
%   %   sampled at 4 MHz
%   load('negtransitionex.mat', 'x', 't');
%   y = falltime(x, 4e6)
%
% See also:  STATELEVELS, RISETIME, SLEWRATE, OVERSHOOT, MIDCROSS.

%   Copyright 2011-2013 The MathWorks, Inc.

% plot if no output specified
plotFlag = nargout==0;

% obtain durations and crossings for all negative-going transitions
[f,lwrCross,uprCross,lwrRef,uprRef] = transdurs(x,-1,plotFlag,varargin{:});

