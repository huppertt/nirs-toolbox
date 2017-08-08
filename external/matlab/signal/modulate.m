function [y,t_out] = modulate(x,Fc,Fs,method,P1)
%MODULATE  Signal modulation for communications simulations.
%   Y = MODULATE(X,Fc,Fs,METHOD,OPT) modulates the message signal X with a
%   carrier frequency Fc and sampling frequency Fs, using the modulation
%   scheme in METHOD.  OPT is an extra sometimes optional parameter whose
%   purpose depends on the modulation scheme you choose.
%
%   Fs must satisfy Fs > 2*Fc + BW, where BW is the bandwidth of the
%   modulated signal.
%
%        METHOD              MODULATION SCHEME
%    'am',      Amplitude modulation, double side-band, suppressed carrier
%    'amdsb-sc' OPT not used.
%    'amdsb-tc' Amplitude modulation, double side-band, transmitted carrier
%               OPT is a scalar which is subtracted from X prior to
%               multiplication by the carrier cosine.  It defaults to
%               min(min(X)) so the offset message signal is positive and
%               has a minimum value of zero.
%    'amssb'    Amplitude modulation, single side-band
%               OPT not used.
%    'fm'       Frequency modulation
%               OPT is a scalar which specifies the constant of frequency
%               modulation kf.  kf = (Fc/Fs)*2*pi/max(max(abs(X))) by
%               default for a maximum frequency excursion of Fc Hertz.
%    'pm'       Phase modulation
%               OPT is a scalar which specifies the constant of phase
%               modulation kp.  kp = pi/max(max(abs(x))) by default for a
%               maximum phase excursion of +/-pi radians.
%    'pwm'      Pulse width modulation
%               If you let OPT = 'centered', the pulses are centered on the
%               carrier period rather than being "left justified".
%    'ppm'      Pulse position modulation
%               OPT is a scalar between 0 and 1 which specifies the pulse
%               width in fractions of the carrier period. It defaults to .1.
%    'qam'      Quadrature amplitude modulation
%               OPT is a matrix the same size as X which is modulated in
%               quadrature with X.
%
%   If X is a matrix, its columns are modulated.
%
%   [Y,T] = MODULATE(...) returns a time vector the same length as Y.
%
%   % Example:
%   %   Frequency modulate a sinusiod with carrier frequency of 3KHz.
%
%   Fs = 8000;                              % sampling frequency
%   t = (0:1000-1)/Fs;                      % time vector
%   s = 4*cos(2*pi*500*t);                  % modulating signal
%   y = modulate(s,3e3,Fs,'fm',0.1);        % modulated signal
%   [ps,fs] = pwelch(s,[],0,length(s),Fs);  % PSD, modulating signal
%   [py,fy] = pwelch(y,[],0,length(y),Fs);  % PSD, modulated signal
%
%   subplot(211);
%   plot(fs/1e3,ps); grid on
%   title('Modulating signal')
%   ylabel('Amplitude(Watts)')
%
%   subplot(212);
%   plot(fy/1e3,py,'r'); grid on
%   title('Modulated Signal')
%   xlabel('Frequency(KHz)');
%   ylabel('Amplitude(Watts)')
%
%   See also DEMOD, VCO, COMM, COMM.MODEM.
% 
%   Note: COMM and COMM.MODEM are in the Communications System Toolbox.

%   Copyright 1988-2014 The MathWorks, Inc.


narginchk(3,5)
% Cast to enforce precision rules
Fc = signal.internal.sigcasttofloat(Fc,'double','modulate','Fc',...
  'allownumeric');
Fs = signal.internal.sigcasttofloat(Fs,'double','modulate','Fs',...
  'allownumeric');
% Checks if X is a valid numeric data input
signal.internal.sigcheckfloattype(x,'','modulate','X');

if Fc >= Fs/2,
    error(message('signal:modulate:AliasedCarrier'));
end

if nargin<4,
    method = 'am';
end

[r,c]=size(x);
if r*c == 0,
    y = []; return
end
if (r==1),   % convert row vector to column
    x = x(:);
    len = c;
else
    len = r;
end

if strcmpi(method,'am') || strcmpi(method,'amdsb-sc'),
    t = (0:1/Fs:((len-1)/Fs))';
    t = t(:,ones(1,size(x,2)));
    y = x.*cos(2*pi*Fc*t);
elseif strcmpi(method,'amdsb-tc'),
    t = (0:1/Fs:((len-1)/Fs))';
    t = t(:,ones(1,size(x,2)));
    if nargin == 5,
      % Cast to enforce precision rules
      offset = signal.internal.sigcasttofloat(P1,'double','modulate',...
        'OPT','allownumeric');
    else
      offset = min(min(x));
    end
    y = ( x - offset ).*cos(2*pi*Fc*t);
elseif strcmpi(method,'amssb'),
    t = (0:1/Fs:((len-1)/Fs))';
    t = t(:,ones(1,size(x,2)));
    y = x.*cos(2*pi*Fc*t) + imag(hilbert(x)).*sin(2*pi*Fc*t);
elseif strcmpi(method,'fm'),
    t = (0:1/Fs:((len-1)/Fs))';
    t = t(:,ones(1,size(x,2)));
    if nargin == 5,
        % Cast to enforce precision rules
        kf = signal.internal.sigcasttofloat(P1,'double','modulate',...
          'OPT','allownumeric');
    else
        x_max = max(max(abs(x)));
        if x_max > 0
            kf = (Fc/Fs)*2*pi/x_max;
            % default- maximum excursion of Fc Hertz
        else
            kf = 0;
        end
    end
    y = cos(2*pi*Fc*t + kf*cumsum(x));   % rectangular integral approx
elseif strcmpi(method,'pm'),
    t = (0:1/Fs:((len-1)/Fs))';
    t = t(:,ones(1,size(x,2)));
    if nargin == 5,
        % Cast to enforce precision rules
        kp = signal.internal.sigcasttofloat(P1,'double','modulate',...
          'OPT','allownumeric');
    else
        kp = pi/max(max(abs(x)));   % default- maximum excursion of +/-pi rads
    end
    y = cos(2*pi*Fc*t + kp*x);
elseif strcmpi(method,'pwm')
    x_max = max(max(x));
    x_min = min(min(x));
    if (x_max > 1) || (x_min < 0)
        error(message('signal:modulate:InvalidRange', 'X', '[0,1]'))
    end
    t = (0:(Fs/Fc*len)-1)';
    if nargin < 5
        P1 = 'left';
    end
    if strcmpi('centered',P1)
        x(len+1,:) = zeros(1,size(x,2));  % catenate a zero at end of input
        % Class(y) must be the same as  class(x) to enforce precision rules
        y = zeros(length(t),size(x,2),class(x));  %#ok<*ZEROLIKE>
        mod_t =  rem(t,(Fs/Fc))/(Fs/Fc);
        ind = floor(t*Fc/Fs+1);
        for i = 1:size(x,2)   % for each column ...
            y(:,i) = (mod_t < x(ind,i)/2)  +  ((1-mod_t) <= x(ind+1,i)/2);
        end
    elseif strcmpi('left',P1)
        % Class(y) must be the same as  class(x) to enforce precision rules
        y = zeros(length(t),size(x,2),class(x));
        for i = 1:size(x,2)   % for each column ...
            y(:,i) = rem(t,(Fs/Fc))/(Fs/Fc) < x(floor(t*Fc/Fs+1),i);
        end
    else
        error(message('signal:modulate:InvalidParam'))
    end
    t = t/Fs;
elseif strcmpi(method,'ptm') || strcmpi(method,'ppm')
    x_max = max(max(x));
    x_min = min(min(x));
    if (x_max > 1) || (x_min < 0)
        error(message('signal:modulate:InvalidRange', 'X', '[0,1]'))
    end
    if nargin<5, P1 = .1; end   % default width
    
    if (P1 >= 1) || (P1 <= 0)
        error(message('signal:modulate:InvalidPDRange'))
    end
    
    if (x_max > 1-P1)
        warning(message('signal:modulate:PulseOverlapped'))
    end
    % Cast to enforce precision rules
    P1 = signal.internal.sigcasttofloat(P1,'double','modulate',...
      'OPT','allownumeric');
    t = (0:(len*Fs/Fc)-1)'/Fs;
    % Class(y) must be the same as  class(x) to enforce precision rules
    y = zeros(length(t)+ceil(P1*Fs/Fc),size(x,2),class(x));
    temp = 0:(P1*Fs/Fc - 1);
    for i=1:size(x,2)    % for each column ...
        ind = ceil( (x(:,i)+(0:len-1)')*Fs/Fc )+1;
        ind = ind(:,ones(1,floor(P1*Fs/Fc) )) + temp(ones(1,length(ind)),:);
        y(ind(:),i) = ones(size(ind(:)));
    end
    y(length(t)+1:length(y),:) = []; % truncate over-long vector
elseif strcmpi(method,'qam')
    if nargin < 5
        error(message('signal:modulate:SignalErr'))
    end
    
    % Checks if X and OPT are valid numeric data inputs
    signal.internal.sigcheckfloattype(x,'','modulate','X');
    signal.internal.sigcheckfloattype(P1,'','modulate','OPT');

    x2 = P1;
    r2 = size(x2, 1);
    if r2 == 1
        x2 = x2.';
    end
    if any(size(x2)~=size(x)),
        error(message('signal:modulate:InvalidDimensions'))
    end
    
    t = (0:1/Fs:((len-1)/Fs))';
    t = t(:,ones(1,size(x,2)));
    y = x.*cos(2*pi*Fc*t) + x2.*sin(2*pi*Fc*t);
end

t = t(:,1);   % only want first column
if (r==1)   % convert y from a column to a row
    y = y.';
    t = t.';
end
if nargout == 2
  if isa(y,'single')
    t_out = single(t);
  else
    t_out = t;
  end
end

