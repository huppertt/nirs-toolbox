function strips(x,sd,Fs,scale)
%STRIPS  Strip plot.
%   STRIPS(X) plots vector X in horizontal strips of length 250.
%   If X is a matrix, STRIPS(X) plots each column of X in horizontal
%   strips.  The left-most column (column 1) is the top horizontal strip.
%
%   STRIPS(X,N) plots vector X in strips that are each N samples long.
%
%   STRIPS(X,SD,Fs) plots vector X in horizontal strips of duration SD
%   seconds given sampling frequency of Fs samples per second.
%
%   STRIPS(X,SD,Fs,SCALE) scales the vertical axes.
%
%   If X is a matrix, STRIPS(X,N), STRIPS(X,SD,Fs), and STRIPS(X,SD,Fs,SCALE)
%   plot the different columns of X on the same strip plot.
%
%   STRIPS ignores the imaginary part of X if it is complex.
%
%   % Example:
%   %   Plot two seconds of a frequency modulated sinusoid in 0.25 
%   %   second strips.
%
%   fs = 1000;                         % Sampling frequency
%   t = 0:1/fs:2;                      % Time vector
%   x = vco(sin(2*pi*t),[10 490],fs);  % FM waveform
%   strips(x,0.25,fs)
%
%   See also PLOT, STEM.

%   Use the syntax strips(X(:),size(X,1),1,SCALE) to get the effect of
%   STRIPS(X) (where X is a matrix) with a SCALE paratmeter.

%   Mark W. Reichelt and Thomas P. Krauss  7-23-93
%   Copyright 1988-2009 The MathWorks, Inc.

error(nargchk(1,4,nargin,'struct'))
if nargin == 1
  if isempty(x), return, end
  Fs = 1;
  scale = 1;
  [sd,c] = size(x);
  if min(sd,c) == 1
      sd = 250;
  end
  x = x(:);
else
  if nargin<3
      Fs = 1;
  end
  if nargin<4
      scale = 1;
  end
  if isempty(x), return, end
  if min(size(x))==1, x = x(:); end   % turn vectors into columns
end

if any(imag(x)~=0), 
  warning(message('signal:strips:ComplexValues'));
  x = real(x);
end

perstrip = ceil(sd * Fs);	% strip duration * number of samples per second
len = size(x,1);                % number of rows
num_sigs = size(x,2);          % number of columns
if rem(len,perstrip) == 1	% leave off last point if it's a singleton
  len = len - 1;
  x = x(1:len,:);
end
num_strips = ceil(len/perstrip);

xmax = max(max( x(~isnan(x)) ));
xmin = min(min( x(~isnan(x)) ));
x0 = 0.5 * (xmin + xmax);

x = scale * x;

% add NaN's to the vector x
NaNind = len+1:perstrip*num_strips;
if ~isempty(NaNind)
    x(NaNind,:)=NaN*ones(length(NaNind),num_sigs);
end

% compute vertical deviation to add to x
del = 0.25 * (xmax-xmin);
sep = (xmax-xmin) + del;
if sep == 0, sep = 1; end
deviation = (num_strips-1:-1:0)*sep;

Y = zeros((perstrip+1)*num_strips,num_sigs);
for l = 1:num_sigs
    y = [reshape(x(:,l),perstrip,num_strips); NaN*ones(1,num_strips)];

    % add vertical deviation to x
    y = y - x0 + deviation(ones(perstrip+1,1),:);
    Y(:,l) = y(:);
end

% compute horizontal (time) axis
t = (0:perstrip-1)'/Fs;
t = t(:,ones(1,num_strips));
t(perstrip+1,:) = NaN + zeros(1,num_strips);
t = t(:);

% compute yticks and yticklabels
yt = (0:num_strips-1)*sep;   % ticks
width = 32;
s = char(ones(num_strips, width) * ' ');
col = width + 1;
for i = 1:num_strips
   str = num2str((i-1)*sd);
   s(i,width-length(str)+1:width) = str;
   col = min(col,width-length(str)+1);
end
s = flipud(s);
s = s(:,col:width);

% plot and set axes properties
newplot;
plot(t,Y)
set(gca,'xlim',[0 sd],'ylim',xmin-x0+[-del sep*num_strips],'ytick',yt,...
    'yticklabel',s,'ygrid','on')

