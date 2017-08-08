function s = square(t,duty)
%SQUARE Square wave generation.
%   SQUARE(T) generates a square wave with period 2*Pi for the
%   elements of time vector T.  SQUARE(T) is like SIN(T), only
%   it creates a square wave with peaks of +1 to -1 instead of
%   a sine wave.
%   SQUARE(T,DUTY) generates a square wave with specified duty
%   cycle. The duty cycle, DUTY, is the percent of the period
%   in which the signal is positive.
%
%   For example, generate a 30 Hz square wave:
%        t = 0:.0001:.0625;
%        y = square(2*pi*30*t);, plot(t,y)
%
%   See also SIN, COS, CHIRP, DIRIC, GAUSPULS, PULSTRAN, RECTPULS,
%   SINC and TRIPULS.

%   Author(s): L. Shure, 2-23-88
%   Copyright 1988-2004 The MathWorks, Inc.

% If no duty specified, make duty cycle 50%.
if nargin < 2
	duty = 50;
end

% Check the input data type. Single precision is not supported.
try
    chkinputdatatype(t,duty);
catch ME
    throwAsCaller(ME);
end

if any(size(duty)~=1),
	error(message('signal:square:SignalErr'))
end

% Compute values of t normalized to (0,2*pi)
tmp = mod(t,2*pi);

% Compute normalized frequency for breaking up the interval (0,2*pi)
w0 = 2*pi*duty/100;

% Assign 1 values to normalized t between (0,w0), 0 elsewhere
nodd = (tmp < w0);

% The actual square wave computation
s = 2*nodd-1;

