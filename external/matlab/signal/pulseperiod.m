function [p,initCross,finalCross,nextCross,midRef] = pulseperiod(x, varargin)
%PULSEPERIOD Period of bilevel waveform pulses
%   P = PULSEPERIOD(X) returns a vector containing the difference between
%   the mid-crossings of the initial transition of each positive-polarity
%   pulse and the next positive-going transition. The signal's sample
%   values are specified by a vector, X, whose sample instants correspond
%   to the index of each sample value. To determine the transitions that
%   compose each pulse, the function estimates the state levels of the
%   input waveform by a histogram method and identifies all regions which
%   cross between a 2% tolerance of each state level.
% 
%   P = PULSEPERIOD(X,Fs) specifies the sample rate, Fs, as a positive
%   scalar, where the first sample instant corresponds to a time of zero.
% 
%   P = PULSEPERIOD(X,T) specifies the sample instants, T, as a vector with
%   the same number of elements as X.
%
%   P = PULSEPERIOD(...,'Polarity',POL) specifies the polarity, POL, of the
%   pulse as either 'positive' or 'negative', where the default value is
%   'positive'. If 'positive' is specified, the function looks for positive
%   polarity pulses (a pulse whose initial transition is positive-going).
%   If 'negative' is specified, the function looks for negative polarity
%   pulses.
% 
%   P = PULSEPERIOD(...,'Tolerance',TOL) specifies the tolerance that the
%   initial and final levels of each transition must be within their
%   respective state levels. The amount, TOL, is a scalar value expressed
%   as percentage of the difference between the upper and lower state
%   levels. The default value is 2.0 (percent).
% 
%   P = PULSEPERIOD(...,'MidPercentReferenceLevel',L) specifies the mid
%   reference level as a percentage of the waveform amplitude, where 0%
%   corresponds to the lower state level, and 100% corresponds to the upper
%   state level. L is a real scalar value. The default value for L is 50
%   (percent).
% 
%   P = PULSEPERIOD(..., 'StateLevels', SL) specifies the levels to use for
%   the lower and upper state levels, SL, as a two-element real row vector
%   whose first and second elements correspond to the lower and upper state
%   levels of the input waveform.
%
%   [P,INITCROSS] = PULSEPERIOD(...) returns a vector, INITCROSS, whose
%   elements correspond to the mid-crossings of the initial transition of
%   each pulse.
%
%   [P,INITCROSS,FINALCROSS] = PULSEPERIOD(...) returns a vector,
%   FINALCROSS, whose elements correspond to the mid-crossings of the final
%   transition of each pulse.
%
%   [P,INITCROSS,FINALCROSS,NEXTCROSS] = PULSEPERIOD(...) returns a vector,
%   NEXTCROSS, whose elements correspond to the mid-crossings of the next
%   detected transition after each pulse.
%
%   [P,INITCROSS,FINALCROSS,NEXTCROSS,MIDLEV] = PULSEPERIOD(...) returns
%   the level, MIDLEV, that corresponds to the mid-reference level.
% 
%   PULSEPERIOD(...) plots the signal and darkens every other identified
%   pulse. It marks the location of the mid crossings, and their associated
%   reference level. The state levels and their associated lower and upper
%   boundaries (adjustable by the TOL parameter) are also plotted.
%
%   % Example 1:
%   %   Plot the pulse period information of a 2.3V pulse waveform 
%   load('pulseex.mat', 'x', 't');
%   pulseperiod(x, t)
%
%   % Example 2:
%   %   Plot the pulse period information of a 2.3V pulse waveform 
%   %   sampled at 4 MHz
%   load('pulseex.mat', 'x', 't');
%   y = pulseperiod(x, 4e6)
%
% See also:  STATELEVELS, MIDCROSS, DUTYCYCLE, PULSEWIDTH, PULSESEP.

%   Copyright 2011-2013 The MathWorks, Inc.

% plot results if no output specified
plotFlag = nargout==0;

% obtain cycle information
[midRef,initCross,finalCross,nextCross] = pulsecycles(x,plotFlag,'PulsePeriod',varargin{:});

% return pulse periods
p = nextCross - initCross;
