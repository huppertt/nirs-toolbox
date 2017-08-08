function [h,ww] = freqs(b,a,w)
%FREQS Laplace-transform (s-domain) frequency response.
%   H = FREQS(B,A,W) returns the complex frequency response vector H
%   of the filter B/A:
%                        nb-1         nb-2
%            B(s)   b(1)s     +  b(2)s     + ... +  b(nb)
%     H(s) = ---- = -------------------------------------
%                        na-1         na-2
%            A(s)   a(1)s     +  a(2)s     + ... +  a(na)
%
%   given the numerator and denominator coefficients in vectors B and A.
%   The frequency response is evaluated at the points specified in
%   vector W (in rad/s).  The magnitude and phase can be graphed by
%   calling FREQS(B,A,W) with no output arguments.
%
%   [H,W] = FREQS(B,A) automatically picks a set of 200 frequencies W on
%   which the frequency response is computed.  FREQS(B,A,N) picks N
%   frequencies.
%
%   % Example 1:
%   %   Find and graph the frequency response of the transfer function
%   %   given by
%   %   H(s) = ( 0.2*s^2 + 0.3*s + 1 )/( s^2 + 0.4s + 1 )
%
%   a = [1 0.4 1];      % Numerator coefficients
%   b = [0.2 0.3 1];    % Denominator coefficients
%   w = logspace(-1,1); % Frequency vector
%   freqs(b,a,w)
%
%   % Example 2:
%   %   Design a fifth-order analog lowpass Bessel filter with an
%   %   approximate constant group delay up to 10,000 rad/s and plot
%   %   the frequency response of the filter using freqs.
%
%   [b,a] = besself(5,10000);   % Bessel analog filter design
%   freqs(b,a)                  % Plot frequency response
%
%   See also LOGSPACE, POLYVAL, INVFREQS, and FREQZ.

% 	Author(s): J.N. Little, 6-26-86
%   	   T. Krauss, 3-19-93, default plots and frequency vector
%   Copyright 1988-2013 The MathWorks, Inc.

narginchk(2,3);
nargoutchk(0,2);

if isempty(b)
  b = 0;
end
if isempty(a)
  a = 0;
end


validateattributes(b,{'numeric'},{'vector'},'freqs','B');
validateattributes(a,{'numeric'},{'vector'},'freqs','A');

[b,a] = removeTrailingZero(b,a);

if nargin == 2,
  w = 200;
end

if length(w) == 1,
  validateattributes(w,{'numeric'},{'scalar','positive','integer'},'freqs','N');
  % Cast to enforce precision rules
  w = double(w);
  wlen = w;
  w_long = freqint(b,a,wlen);
  % need to interpolate long frequency vector:
  xi = linspace(1,length(w_long),wlen).';
  w = 10.^interp1(1:length(w_long),log10(w_long),xi,'linear');
else
  validateattributes(w,{'numeric'},{'vector'},'freqs','W');
  % Cast to enforce precision rules
  w = double(w);
end

s = 1i*w;
hh = polyval(b,s) ./ polyval(a,s);

if nargout == 0,
  % make sure W is never decreasing
  if w(1) > w(length(w))
    error(message('signal:freqs:IncreasingW','W','W'));
  end
  newplot;
  mag = abs(hh);   phase = angle(hh)*180/pi;
  subplot(211),loglog(w,mag),set(gca,'xgrid','on','ygrid','on')
  set(gca,'xlim',[w(1) w(length(w))])
  xlabel(getString(message('signal:freqs:Frequencyrads')))
  ylabel(getString(message('signal:freqs:Magnitude')))
  ax = gca;
  subplot(212), semilogx(w,phase),set(gca,'xgrid','on','ygrid','on')
  set(gca,'xlim',[w(1) w(length(w))])
  xlabel(getString(message('signal:freqs:Frequencyrads')))
  ylabel(getString(message('signal:freqs:Phasedegrees')))
  axes(ax) 
elseif nargout == 1,
  h = hh;
elseif nargout == 2,
  h = hh; 
  % Cast to enforce precision rules
  if isa(h,'single')
    ww = single(w);
  else
    ww = w;
  end
end
% end freqs

function w=freqint(a,b,npts)
%FREQINT Auto-ranging algorithm for Bode frequency response
%   W=FREQINT(NUM,DEN,Npts)

% Generate more points where graph is changing rapidly.
% Calculate points based on eigenvalues and transmission zeros.

ep=roots(b);
tz=roots(a);

if isempty(ep), ep=-1000; end

% Note: this algorithm does not handle zeros greater than 1e5
ez=[ep(imag(ep)>=0);tz(abs(tz)<1e5&imag(tz)>=0)];

% Round first and last frequencies to nearest decade
integ = abs(ez)<1e-10; % Cater for systems with pure integrators
highfreq=round(log10(max(3*abs(real(ez)+integ)+1.5*imag(ez)))+0.5);
lowfreq=round(log10(0.1*min(abs(real(ez+integ))+2*imag(ez)))-0.5);

% Define a base range of frequencies
diffzp=length(ep)-length(tz);
w=logspace(lowfreq,highfreq,npts+diffzp+10*(sum(abs(imag(tz))<abs(real(tz)))>0));
ez=ez(imag(ez)>abs(real(ez)));

% Oscillatory poles and zeros
if ~isempty(ez)
  f=w;
  npts2=2+8/ceil(abs((diffzp+eps)/10));
  [dum,ind]=sort(-abs(real(ez))); %#ok
  z=[];
  for i=ind'
    r1=max([0.8*imag(ez(i))-3*abs(real(ez(i))),10^lowfreq]);
    r2=1.2*imag(ez(i))+4*abs(real(ez(i)));
    z=z(z>r2|z<r1);
    f=f(f>r2|f<r1);
    z=[z,logspace(log10(r1),log10(r2),sum(w<=r2&w>=r1)+npts2)]; %#ok<AGROW>
  end
  w=sort([f,z]);
end
w = w(:);
% Cast to enforce precision rules
w = double(w);

% end freqint

function [b,a] = removeTrailingZero(b,a)

b_len = numel(b);
a_len = numel(a);
b_lastnzidx = find(b,1,'last');
a_lastnzidx = find(a,1,'last');
trz_len = min(b_len-b_lastnzidx,a_len-a_lastnzidx);
if trz_len > 0
  b = b(1:b_len-trz_len);
  a = a(1:a_len-trz_len);
end


% end removeTrailingZero

