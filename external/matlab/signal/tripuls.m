function y = tripuls(t,Tw,skew)
%TRIPULS Sampled aperiodic triangle generator.
%   TRIPULS(T) generates samples of a continuous, aperiodic,
%   unity-height triangle at the points specified in array T, centered
%   about T=0.  By default, the triangle is symmetric and has width 1.
%
%   TRIPULS(T,W) generates a triangle of width W.
%
%   TRIPULS(T,W,S) allows the triangle skew S to be adjusted.  The
%   skew parameter must be in the range -1 < S < +1, where 0 generates
%   a symmetric triangle.
%
%   % Example:
%   %   Create a triangular pulse with width 0.4.
%
%   fs = 10000;         % Sampling frequency (samples/sec)
%   t = -1:1/fs:1;      % Time Vector
%   w = .4;             % Triangle Width
%   x = tripuls(t,w);   % Sampled aperiodic triangle
%   plot(t,x); 
%   xlabel('Time (sec)');ylabel('Amplitude');
% 
%   See also SQUARE, SIN, COS, CHIRP, DIRIC, GAUSPULS, PULSTRAN and
%   RECTPULS.

%   Copyright 1988-2013 The MathWorks, Inc.

narginchk(1,3);

if nargin<2, Tw=1;   end
if nargin<3, skew=0; end
if (skew<-1) || (skew>1),
  error(message('signal:tripuls:InvalidRange'));
end

% Cast to enforce precision rules
t = signal.internal.sigcasttofloat(t,'double','tripuls','T','allownumeric');
Tw = signal.internal.sigcasttofloat(Tw,'double','tripuls','W',...
  'allownumeric');
skew = signal.internal.sigcasttofloat(skew,'double','tripuls','S',...
  'allownumeric');

% Compute triangle function output:
Ta=Tw/2*skew;
y=zeros(size(t));
i = find( (t>(-Tw/2)) & (t<Ta) );
y(i) = (2*t(i)+Tw)./(Tw*(1+skew));
i = find( (t>Ta) & (t<Tw/2) );
y(i) = 1 - (2*t(i)-skew*Tw)./(Tw*(1-skew));
y(t==Ta) = 1.0;

