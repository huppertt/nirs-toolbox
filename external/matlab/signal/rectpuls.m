function y = rectpuls(t,Tw)
%RECTPULS Sampled aperiodic rectangle generator.
%   RECTPULS(T) generates samples of a continuous, aperiodic,
%   unity-height rectangle at the points specified in array T, centered
%   about T=0.  By default, the rectangle has width 1.  Note that the
%   interval of non-zero amplitude is defined to be open on the right,
%   i.e., RECTPULS(-0.5)=1 while RECTPULS(0.5)=0.
%
%   RECTPULS(T,W) generates a rectangle of width W.
%
%   % Example:
%   %   Generate 2 seconds of a rectangular pulse with a sample rate of 
%   %   10Khz and a width of 20ms.
%
%   fs = 10000;                             % Number of samples per second
%   t = -1:1/fs:1;                          % Time vector
%   x2 = rectpuls(t,20e-3);                 % Generating rectangular pulse 
%   plot(t,x2), axis([-0.1 0.1 -0.2 1.2])   % Plotting the pulse
%
%   See also GAUSPULS, TRIPULS, PULSTRAN.

%   Copyright 1988-2013 The MathWorks, Inc.

narginchk(1,2);
if nargin<2, Tw=1;   end

% Cast to enforce precision rules
t = signal.internal.sigcasttofloat(t,'double','rectpuls','T',...
  'allownumeric');

Tw = signal.internal.sigcasttofloat(Tw,'double','rectpuls','W',...
  'allownumeric');

% Returns unity in interval [-Tw/2,+Tw/2) (right side of interval is open)
% Because of numerical concerns, we use eps as the tolerance
y = double(abs(t)<Tw/2-eps);
y(abs(t-(-Tw/2))<eps) = 1.0;

% end of rectpuls.m
