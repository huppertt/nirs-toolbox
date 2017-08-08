function [w,initCross,finalCross,midRef] = pulsewidth(x, varargin)
%PULSEWIDTH Width of bilevel waveform pulses
%   W = PULSEWIDTH(X) returns a vector containing the difference between
%   the mid-crossings of the initial and final transitions of each
%   positive-polarity pulse found in the input signal. The signal's sample
%   values are specified by a vector, X, whose sample instants correspond
%   to the index of each sample value. To determine the transitions, the
%   function estimates the state levels of the input waveform by a
%   histogram method and identifies all regions which cross between a 2%
%   tolerance of each state level.
% 
%   W = PULSEWIDTH(X,Fs) specifies the sample rate, Fs, as a positive
%   scalar, where the first sample instant corresponds to a time of zero.
% 
%   W = PULSEWIDTH(X,T) specifies the sample instants, T, as a vector with
%   the same number of elements as X.
%
%   W = PULSEWIDTH(...,'Polarity',POL) specifies the polarity, POL, of the
%   pulse as either 'positive' | 'negative', where the default value is
%   'positive'. If 'positive' is specified, the function looks for positive
%   polarity pulses (a pulse whose initial transition is positive-going).
%   If 'negative' is specified, the function looks for negative polarity
%   pulses.
% 
%   W = PULSEWIDTH(...,'Tolerance',TOL) specifies the tolerance that the
%   initial and final levels of each transition must be within their
%   respective state levels. The amount, TOL, is a real scalar expressed as
%   the percentage of the difference between the upper and lower state
%   levels. The default value is 2.0 (percent).
% 
%   W = PULSEWIDTH(...,'MidPercentReferenceLevel',L) specifies the mid
%   reference level as a percentage of the waveform amplitude, where 0%
%   corresponds to the lower state level, and 100% corresponds to the upper
%   state level. L is a real scalar. The default value for L is 50
%   (percent).
% 
%   W = PULSEWIDTH(...,'StateLevels',SL) specifies the levels to use for
%   the lower and upper state levels, SL, as a two-element real row vector
%   whose first and second elements correspond to the lower and upper state
%   levels of the input waveform.
% 
%   [W,INITCROSS] = PULSEWIDTH(...) returns a vector, INITCROSS, whose
%   elements correspond to the mid-crossings of the initial transition of
%   each pulse.
% 
%   [W,INITCROSS,FINALCROSS] = PULSEWIDTH(...) returns a vector, FINALCROSS,
%   whose elements correspond to the mid-crossings of the final transition 
%   of each pulse.
%
%   [W,INITCROSS,FINALCROSS,MIDLEV] = PULSEWIDTH(...) returns the
%   level, MIDLEV, that corresponds to the mid-reference level.
%
%   PULSEWIDTH(...) plots the signal and darkens the regions of each pulse
%   where pulse width is computed. It marks the location of the mid
%   crossings, and their associated reference level. The state levels and
%   their associated lower and upper boundaries (adjustable by the TOL
%   parameter) are also plotted.
%
%   % Example 1:
%   %   Compute the pulse width information of a 2.3V pulse waveform 
%   load('pulseex.mat', 'x', 't');
%   pulsewidth(x, t)
%
%   % Example 2:
%   %   Compute the pulse width information of a 2.3V pulse waveform 
%   %   sampled at 4 MHz
%   load('pulseex.mat', 'x', 't');
%   y = pulsewidth(x, 4e6)
%
% See also:  STATELEVELS, MIDCROSS, DUTYCYCLE, PULSESEP, PULSEPERIOD.

%   Copyright 2011-2013 The MathWorks, Inc.

% plot results if no output specified
plotFlag = nargout==0;

% obtain cycle information (ignoring edge of next pulse)
[midRef,initCross,finalCross] = pulsecycles(x,plotFlag,'PulseWidth',varargin{:});

% return pulse width
w = finalCross - initCross;
