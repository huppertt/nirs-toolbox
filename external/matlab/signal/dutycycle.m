function [d,initCross,finalCross,nextCross,midRef] = dutycycle(x, varargin)
%DUTYCYCLE Duty cycle of bilevel waveform pulses
%   D = DUTYCYCLE(X) returns the ratio of the pulse width to the pulse
%   period for each positive-polarity pulse. D is a vector whose number of
%   elements correspond to the number of pulse periods found in the input
%   signal. The signal's sample values are specified by a vector, X, whose
%   sample instants correspond to the index of each sample value. If no
%   complete pulse periods are found, D is an empty vector. To determine
%   the transitions that compose each pulse, the function estimates the
%   state levels of the input waveform by a histogram method and identifies
%   all regions which cross between a 2% tolerance of each state level.
% 
%   D = DUTYCYCLE(X,Fs) specifies the sample rate, Fs, as a positive
%   scalar, where the first sample instant corresponds to a time of zero.
%
%   D = DUTYCYCLE(X,T) specifies the sample instants, T, as a vector with
%   the same number of elements as X.
% 
%   D = DUTYCYCLE(TAU,PRF) returns the duty cycle D for given real scalar
%   input pulse width TAU (in seconds) and pulse repetition frequency PRF
%   (in Hz). The product of TAU and PRF must be less than or equal to 1.
%
%   D = DUTYCYCLE(...,'Polarity',POL) specifies the polarity, POL, of the
%   pulse as either 'positive' or 'negative', where the default value is
%   'positive'. If 'positive' is specified, then the function looks for
%   positive polarity pulses (a pulse whose initial transition is
%   positive-going). If 'negative' is specified, then the function looks
%   for negative polarity pulses.
% 
%   D = DUTYCYCLE(...,'Tolerance',TOL) specifies the tolerance that the
%   initial and final levels of each transition must be within their
%   respective state levels. The amount, TOL, is a scalar value expressed
%   as percentage of the difference between the upper and lower state
%   levels. The default value is 2.0 (percent)
% 
%   D = DUTYCYCLE(...,'MidPercentReferenceLevel',L) specifies the mid
%   reference level as a percentage of the waveform amplitude, where 0%
%   corresponds to the lower state level, and 100% corresponds to the upper
%   state level. Specify L as a scalar value. The default value for L is 50
%   (percent).
% 
%   D = DUTYCYCLE(...,'StateLevels',SL) specifies the levels to use for the
%   lower and upper state levels, SL, as a two-element real row vector
%   whose first and second element correspond to the lower and upper state
%   levels of the input waveform.
% 
%   [D,INITCROSS] = DUTYCYCLE(...) returns a vector, INITCROSS, whose
%   elements correspond to the mid-crossings of the initial transition of
%   each pulse.
%
%   [D,INITCROSS,FINALCROSS] = DUTYCYCLE(...) returns a vector, FINALCROSS,
%   whose elements correspond to the mid-crossings of the final transition
%   of each pulse.
%
%   [D,INITCROSS,FINALCROSS,NEXTCROSS] = DUTYCYCLE(...) returns a vector,
%   NEXTCROSS, whose elements correspond to the mid-crossings of the next
%   detected transition after each pulse.
%
%   [D,INITCROSS,FINALCROSS,NEXTCROSS,MIDLEV] = DUTYCYCLE(...) returns the
%   scalar, MIDLEV, that corresponds to the mid-reference level.
% 
%   DUTYCYCLE(...) plots the signal and marks the location of the mid
%   crossings, and their associated reference level. The state levels and
%   their associated lower and upper boundaries (adjustable by the TOL
%   parameter) are also plotted.
%
%   % Example 1:
%   %   Compute the duty cycle information of a 2.3V step waveform
%   load('pulseex.mat', 'x', 't');
%   y = dutycycle(x, t)
%
%   % Example 2:
%   %   Compute the duty cycle information of a 2.3V step waveform 
%   %   sampled at 4 MHz
%   load('pulseex.mat', 'x', 't');
%   y = dutycycle(x, 4e6)
%
% See also:  STATELEVELS, MIDCROSS, PULSEWIDTH, PULSESEP, PULSEPERIOD.

%   Copyright 2011-2013 The MathWorks, Inc.

%#ok<*EMCLS>
%#ok<*EMCA>
%#codegen

if isscalar(x)
    cond = nargin==2;
    if ~cond
        coder.internal.assert(cond,'signal:dutycycle:TauPrfWrongNumberOfInputs');
    end
    cond = nargout<=1;
    if ~cond
        coder.internal.assert(cond,'signal:dutycycle:TauPrfTooManyOutputs');
    end
    prf = varargin{1};
    validateattributes(x,{'double'},{'scalar','positive','finite'},'dutycycle','TAU');
    validateattributes(prf,{'double'},{'positive','finite'},'dutycycle','PRF');
    d = x*prf;

    cond = any(d>1);
    if cond
        coder.internal.errorIf(cond,'signal:dutycycle:TauPrfExceededUnity');
    end
  
else

    if ~coder.target('MATLAB')
    end
    
    % plot results if no output specified
    plotFlag = nargout==0;

    % obtain cycle information
    [midRef,initCross,finalCross,nextCross] = pulsecycles(x,plotFlag,'DutyCycle',varargin{:});

    % return pulsewidth / pulseperiod
    d = (finalCross - initCross) ./ (nextCross - initCross);
end
 

