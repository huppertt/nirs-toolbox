function y=pulstran(t,d,func,varargin)
%PULSTRAN Pulse train generator.
%   PULSTRAN generates pulse trains from either continuous functions
%   or sampled prototype pulses.
%
%   CONTINUOUS FUNCTIONS
%   Y=PULSTRAN(T,D,'func') generates a pulse train based on samples of
%   a continuous function, 'func'.  The function is evaluated over
%   the range of argument values specified in array T, after removing
%   a scalar argument offset taken from the vector D.  Thus, the
%   function is evaluated length(D) times, and the sum of the
%   evaluations Y = func(t-D(1)) + func(t-D(2)) + ... is returned.
%   Note that 'func' must be a vectorized function which can take an
%   array T as an argument.
%
%   Y=PULSTRAN(T,D,@func) may be used instead of PULSTRAN(T,D,'func')
%   to specify a function handle instead of a string function name.
%
%   An optional gain factor may be applied to each delayed evaluation
%   by specifying D as a 2-column matrix, with the offset defined in
%   column 1 and associated gain in column 2 of D.  Note that a row
%   vector will be interpreted as specifying argument delays only.
%
%   PULSTRAN(T,D,'func',P1,P2,...) allows additional parameters to be
%   passed to 'func' as necessary, e.g., the function is called as
%   func( T-D(1), P1, P2, ... ), etc.
%
%   SAMPLED PROTOTYPE PULSES
%   PULSTRAN(T,D,P,FS) generates a pulse train consisting of the sum
%   of multiple delayed interpolations of the prototype pulse
%   specified in vector P sampled at the rate FS.  T and D are defined
%   as above.  It is assumed that P spans the time interval [0,
%   (length(P)-1)/FS], and its samples are identically zero outside
%   this interval.  By default, linear interpolation is used for
%   generating delays.
%
%   PULSTRAN(T,D,P) assumes that FS=1 Hz, and PULSTRAN(..., 'method')
%   specifies alternative interpolation methods.  See INTERP1 for a
%   list of available methods.
%
%   EXAMPLES
%
%   % Example 1: Generate an asymmetric sawtooth waveform with a
%                % repetition frequency of 3 Hz and a sawtooth width of 0.1 sec.
%                % The signal is to be 1 second long with a sample rate of 1kHz.
%
%                T = 0 : 1/1E3 : 1;  % 1 kHz sample freq for 1 sec
%                D = 0 : 1/3 : 1;    % 3 Hz repetition freq
%                Y = pulstran(T,D,'tripuls',0.1,-1); plot(T,Y)
%
%   % Example 2: Generate a periodic Gaussian pulse signal at 10 kHz,
%                % with 50% bandwidth.  The pulse repetition frequency is 1 kHz,
%                % sample rate is 50 kHz, and pulse train length is 10msec.  The
%                % repetition amplitude should attenuate by 0.8 each time.  Uses
%                % a function handle to refer to the generator function.
%
%                T = 0 : 1/50E3 : 10E-3;
%                D = [0 : 1/1E3 : 10E-3 ; 0.8.^(0:10)]';
%                Y = pulstran(T,D,@gauspuls,10E3,.5); plot(T,Y)
%
%
%   See also GAUSPULS, RECTPULS, TRIPULS, SINC.

%   Author(s): D. Orofino, 4/96
%   Copyright 1988-2008 The MathWorks, Inc.

if nargin<3,
  error(message('signal:pulstran:Nargchk'));
end
% Check that d is either a vector, or an Mx2 matrix.
[dm,dn]=size(d);
if dm>1 && dn>2,
  error(message('signal:pulstran:InvalidDimensionsD', 'D', 'Mx2'));
end
% Force d to a column vector if it is a vector.  This assumes that
% d=[q1 q2] implies 2 delays and not one delay with an amplitude.
if dm==1, d=d.'; end

% Check the input data type. Single precision is not supported.
try
    chkinputdatatype(t,d);
catch ME
    throwAsCaller(ME);
end
  
y=zeros(size(t));  % Allocate result

if ischar(func) || isa(func,'function_handle'),
  % Continuous function (specified as a string function name
  % or as a function handle):
  
  if size(d,2)==1,
    % Unity amplitudes:
    for i=1:length(d),
      y = y + feval(func,t-d(i),varargin{:});
    end
  else
    % Arbitrary amplitudes:
    for i=1:length(d),
      y = y + d(i,2).*feval(func,t-d(i,1),varargin{:});
    end
  end

else
  % Sampled prototype:
  if nargin==5,
    % (...,fs,'method') specified:
    fs = varargin{1};
    method = varargin{2};

  elseif nargin==4,
    % (...,fs) or (...,'method') specified
    method = varargin{1};
    if ischar(method),
      % method specified:
      fs = 1;
    else
      % fs specified:
      fs = method;
      method = 'linear';
    end
  else
    % Default both fs and 'method':
    method = 'linear';
    fs = 1;
  end

  % Check the input data type. Single precision is not supported.
  try
      chkinputdatatype(fs);
  catch ME
      throwAsCaller(ME);
  end

  P = func;                % Get prototype pulse
  if ~isvector(P),
    error(message('signal:pulstran:InvalidDimensionsP', 'P'));
  end
  x = (-3:length(P)+2)/fs;  % Source sample time instants
  P = [zeros(3,1);P(:);zeros(3,1)];
  xi = t(:);                % Must be vector for interp1
  tsiz = size(t);           % Original shape of t array

  if size(d,2)==1,
    % Unity amplitudes:
    for i=1:length(d),
      Dy = interp1(x,P,xi-d(i,1),method);
      Dy(find(isnan(Dy))) = 0;
      y = y + reshape(Dy,tsiz);
    end
  else
    % Arbitrary amplitudes:
    for i=1:length(d),
      Dy = interp1(x,P,xi-d(i,1),method);
      Dy(find(isnan(Dy))) = 0;
      y = y + d(i,2).*reshape(Dy,tsiz);
    end
  end

end

% end of pulstran.m

function y=isvector(x,dim)
%ISVECTOR True for a vector in N-D space.  ISVECTOR(X) returns 1 is X
%   has at most 1 dimension with size greater than 1. 
%
%   ISVECTOR(X,DIM) is true if X is a vector along the dimension DIM.
%
%   See also ISREAL.

sx = size(x)~=1;   % Dimensions w/size~=0
y = (sum(sx)<=1);  % Is it a vector?
if nargin>1 && y,
  y = find(sx)==dim;
end

% end of ISVECTOR.M
