function [s,initCross,finalCross,nextCross,midRef] = pulsesep(x, varargin)
%PULSESEP Separation between bilevel waveform pulses
%   S = PULSESEP(X) returns a vector containing the differences between the
%   mid-crossings of each final negative-going transition of every
%   positive-polarity pulse and the next positive-going transition. The
%   signal's sample values are specified by a vector, X, whose sample
%   instants correspond to the index of each sample value. To determine the
%   transitions that compose each pulse, the function estimates the state
%   levels of the input waveform by a histogram method and identifies all
%   regions which cross between a 2% tolerance of each state level.
% 
%   S = PULSESEP(X,Fs) specifies the sample rate, Fs, as a positive
%   scalar, where the first sample instant corresponds to a time of zero.
% 
%   S = PULSESEP(X,T) specifies the sample instants, T, as a vector with
%   the same number of elements as X.
%
%   S = PULSESEP(...,'Polarity',POL) specifies the polarity, POL, of the
%   pulse as either 'positive' or 'negative', where the default value is
%   'positive'. If 'positive' is specified, the function looks for positive
%   polarity pulses (a pulse whose initial transition is positive-going).
%   If 'negative' is specified, the function looks for negative polarity
%   pulses.
% 
%   S = PULSESEP(...,'Tolerance', TOL) specifies the tolerance that the
%   initial and final levels of each transition must be within their
%   respective state levels. The amount, TOL, is a scalar value expressed
%   as percentage of the difference between the upper and lower state
%   levels. The default value is 2.0 (percent).
% 
%   S = PULSESEP(...,'MidPercentReferenceLevel',L) specifies the mid
%   reference level as a percentage of the waveform amplitude, where 0%
%   corresponds to the lower state level, and 100% corresponds to the upper
%   state level. L is a real scalar value. The default value for L is 50
%   (percent).
% 
%   S = PULSESEP(...,'StateLevels',SL) specifies the levels to use for the
%   lower and upper state levels, SL, as a two-element real row vector
%   whose first and second elements correspond to the lower and upper state
%   levels of the input waveform.
% 
%   [S,INITCROSS] = PULSESEP(...) returns a vector, INITCROSS, whose
%   elements correspond to the mid-crossings of the initial transition of
%   each pulse.
%
%   [S,INITCROSS,FINALCROSS] = PULSESEP(...) returns a vector, FINALCROSS,
%   whose elements correspond to the mid-crossings of the final transition
%   of each pulse.
%
%   [S,INITCROSS,FINALCROSS,NEXTCROSS] = PULSESEP(...) returns a vector,
%   NEXTCROSS, whose elements correspond to the mid-crossings of the next
%   detected transition after each pulse.
%
%   [S,INITCROSS,FINALCROSS,NEXTCROSS,MIDLEV] = PULSESEP(...) returns the
%   level, MIDLEV, that corresponds to the mid-reference level.
%
%   PULSESEP(...) plots the signal and darkens the regions between each
%   pulse where pulse separation is computed. It marks the location of the
%   mid crossings, and their associated reference level. The state levels
%   and their associated lower and upper boundaries (adjustable by the TOL
%   parameter) are also plotted.
%
%   % Example 1:
%   %   Compute the pulse separation information of a 2.3V pulse waveform 
%   load('pulseex.mat', 'x', 't');
%   pulsesep(x, t)
%
%   % Example 2:
%   %   Compute the pulse separation information of a 2.3V pulse waveform 
%   %   sampled at 4 MHz
%   load('pulseex.mat', 'x', 't');
%   y = pulsesep(x, 4e6)
%
% See also:  STATELEVELS, MIDCROSS, DUTYCYCLE, PULSEWIDTH, PULSEPERIOD.

%   Copyright 2011-2013 The MathWorks, Inc.

% plot results if no output specified
plotFlag = nargout==0;

% obtain cycle information
[midRef,initCross,finalCross,nextCross] = pulsecycles(x,plotFlag,'PulseSeparation',varargin{:});

% return pulse separation
s = nextCross - finalCross;
