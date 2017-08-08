function [s,lwrCross,uprCross,lwrRef,uprRef] = slewrate(x, varargin)
%SLEWRATE Slew rate of bilevel waveform transitions
%   S = SLEWRATE(X) returns a vector of the ratios of the level difference
%   to the time duration between the points where each transition crosses
%   the 10% and 90% reference levels. The signal's sample values are
%   specified by a vector, X, whose sample instants correspond to the index
%   of each sample value. To determine the transitions, the function
%   estimates the state levels of the input waveform by a histogram method
%   and identifies all regions which cross between a 2% tolerance of each
%   state level.
% 
%   S = SLEWRATE(X,Fs) specifies the sample rate, Fs, as a positive scalar,
%   where the first sample instant corresponds to a time of zero.  
%
%   S = SLEWRATE(X,T) specifies the sample instants, T, as a vector with
%   the same number of elements as X.
%
%   S = SLEWRATE(...,'Tolerance',TOL) specifies the tolerance that the
%   initial and final levels of each transition must be within their
%   respective state levels. The amount, TOL, is a scalar expressed
%   as the percentage of the difference between the upper and lower state
%   levels. The default value is 2.0 (percent).
%
%   S = SLEWRATE(...,'PercentReferenceLevels',PRL) specifies the reference
%   levels as a percentage of the waveform amplitude, where 0% corresponds
%   to the lower state level, and 100% corresponds to the upper state
%   level. PRL is a two-element real row vector whose elements correspond
%   to the lower and upper percent reference levels. The default value for
%   this parameter is [10 90].
%
%   S = SLEWRATE(...,'StateLevels',SL) specifies the levels to use for the
%   lower and upper state levels, SL, as a two-element real row vector
%   whose first and second elements correspond to the lower and upper state
%   levels of the input waveform.
%
%   [S,LT,UT] = SLEWRATE(...) returns vectors, LT and UT, whose elements
%   correspond to the time instants where X crosses the lower and upper
%   percent reference levels.
%
%   [S,LT,UT,LL,UL] = SLEWRATE(...) returns the levels, LL and UL, that
%   correspond to the lower and upper percent reference levels.
%
%   SLEWRATE(...) plots the signal and darkens the regions of each
%   transition where slew rate is computed. It marks the lower and upper
%   crossings, and their associated reference levels. The state levels and
%   their associated lower and upper boundaries (adjustable by the TOL
%   parameter) are also plotted.
%
%   % Example 1:
%   %   Plot the slew rate information of a 2.3V step waveform 
%   load('transitionex.mat', 'x', 't');
%   slewrate(x, t)
%
%   % Example 2:
%   %   Compute the slew rate information of a 2.3V step waveform 
%   %   sampled at 4 MHz
%   load('transitionex.mat', 'x', 't');
%   y = slewrate(x, 4e6)
%
% See also:  STATELEVELS, RISETIME, FALLTIME, OVERSHOOT, MIDCROSS.

%   Copyright 2011-2013 The MathWorks, Inc.

% plot if no output specified
plotFlag = nargout==0;

% obtain crossings for all transitions
[~,lwrCross,uprCross,lwrRef,uprRef] = transdurs(x,0,plotFlag,varargin{:});

% return slew rate 
s = (uprRef - lwrRef) ./ (uprCross - lwrCross);

