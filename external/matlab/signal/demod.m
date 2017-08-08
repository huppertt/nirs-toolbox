function [x,x2] = demod(y, Fc, Fs, method, P1)
%DEMOD Signal demodulation for communications simulations.
%   X = DEMOD(Y,Fc,Fs,METHOD,OPT) demodulates the carrier signal Y with a 
%   carrier frequency Fc and sampling frequency Fs, using the demodulation
%   scheme in METHOD.  OPT is an extra, sometimes optional, parameter whose
%   purpose depends on the demodulation scheme you choose.
%
%   Fs must satisfy Fs > 2*Fc + BW, where BW is the bandwidth of the
%   modulated signal.
%
%        METHOD            DEMODULATION SCHEME
%    'am',      Amplitude demodulation, double side-band, suppressed carrier
%    'amdsb-sc' OPT not used.
%    'amdsb-tc' Amplitude demodulation, double side-band,transmitted carrier
%               OPT is a scalar which is subtracted from the decoded message
%               signal.  It defaults to zero.
%    'amssb'    Amplitude demodulation, single side-band
%               OPT not used.
%    'fm'       Frequency demodulation
%               OPT is a scalar which specifies the constant of frequency 
%               modulation kf, which defaults to 1.
%    'pm'       Phase demodulation
%               OPT is a scalar which specifies the constant of phase 
%               modulation kp, which defaults to 1.
%    'pwm'      Pulse width demodulation
%               By setting OPT = 'centered' you tell DEMOD that the pulses 
%               are centered on the carrier period rather than being 
%               "left justified".
%    'ppm'      Pulse position demodulation
%               OPT is not used.
%    'qam'      Quadrature amplitude demodulation
%               For QAM signals, use [X1,X2] = DEMOD(Y,Fc,Fs,'qam')
%
%   If Y is a matrix, its columns are demodulated.
%
%   % Example:
%   %   Demodulate a received signal which is frequency modulated at 
%   %   carrier frequency 3KHz.
%   
%   Fs = 8000;                              % sampling frequency
%   t = (0:1000-1)/Fs;                      % time vector
%   s = 4*cos(2*pi*500*t);                  % modulating signal
%   x = modulate(s,3e3,Fs,'fm',0.1);        % modulated signal
%   rx = x + sqrt(1e-3)*randn(size(x));     % received signal
%   y = demod(rx,3e3,Fs,'fm');              % demodulated signal
%   [px,fx] = pwelch(rx,[],0,length(x),Fs); % PSD, received signal
%   [py,fy] = pwelch(y,[],0,length(y),Fs);  % PSD, demodulated signal
% 
%   subplot(211);
%   plot(fx/1e3,px); grid on
%   title('Received signal')
%   ylabel('Amplitude(Watts)')
% 
%   subplot(212);
%   plot(fy/1e3,py,'r'); grid on
%   title('Demodulated Signal')
%   xlabel('Frequency(KHz)');
%   ylabel('Amplitude(Watts)')
%
%   See also MODULATE PAMDEMOD QAMDEMOD GENQAMDEMOD FSKDEMOD
%            PSKDEMOD MSKDEMOD.
%
%   Note:  PAMDEMOD, QAMDEMOD, GENQAMDEMOD, FSKDEMOD, PSKDEMOD and MSKDEMOD
%          are in the Communications System Toolbox.

%   Copyright 1988-2014 The MathWorks, Inc.

narginchk(3,5)
% Cast to enforce precision rules
Fc = signal.internal.sigcasttofloat(Fc,'double','demod','Fc',...
  'allownumeric');
Fs = signal.internal.sigcasttofloat(Fs,'double','demod','Fs',...
  'allownumeric');

% Checks if Y is a valid numeric data input. Cast to double since filtfilt
% does not accept single inputs. 
isInputSingle = signal.internal.sigcheckfloattype(y,'single','demod','Y');
y = double(y);

if Fc >= Fs/2,
	error(message('signal:demod:InvalidRange'));
end

if nargin<4,
    method = 'am';
end

[r,c]=size(y);
if r*c == 0,
    x = []; return
end
if (r==1),   % convert row vector to column
    y = y(:);  len = c;
else
    len = r;
end

if strcmpi(method,'am') || strcmpi(method,'amdsb-sc') || ...
        strcmpi(method,'amdsb-tc') || strcmpi(method,'amssb'),
    t = (0:1/Fs:((len-1)/Fs))';
    t = t(:,ones(1,size(y,2)));
    x = y.*cos(2*pi*Fc*t);
    [b,a]=butter(5,Fc*2/Fs);
    for i = 1:size(y,2),
        x(:,i) = filtfilt(b,a,x(:,i));
    end
    if strcmpi(method,'amdsb-tc'),
        if nargin < 5,
            P1 = 0;
        end
        P1 = signal.internal.sigcasttofloat(P1,'double','demod','OPT',...
          'allownumeric');        
        x = x - P1;
    end
%    if strcmpi(method,'amdsb-sc')|strcmpi(method,'am'),
%        x = x/2;
%    end
elseif strcmpi(method,'fm'),
    if nargin < 5, P1 = 1; end
    P1 = signal.internal.sigcasttofloat(P1,'double','demod','OPT',...
      'allownumeric');
    t = (0:1/Fs:((len-1)/Fs))';
    t = t(:,ones(1,size(y,2)));
    yq = hilbert(y).*exp(-1i*2*pi*Fc*t);
    x = (1/P1)*[zeros(1,size(yq,2)); diff(unwrap(angle(yq)))];
elseif strcmpi(method,'pm'),
    if nargin < 5, P1 = 1; end
    P1 = signal.internal.sigcasttofloat(P1,'double','demod','OPT',...
      'allownumeric');
    t = (0:1/Fs:((len-1)/Fs))';
    t = t(:,ones(1,size(y,2)));
    yq = hilbert(y).*exp(-1i*2*pi*Fc*t);
    x = (1/P1)*angle(yq);
elseif strcmpi(method,'pwm'),
    % precondition input by thresholding:
    y = y>.5;
    t = (0:1/Fs:((len-1)/Fs))';
    len = ceil( len * Fc / Fs);   % length of message signal
    x = zeros(len,size(y,2));
    if nargin < 5
        P1 = 'left';
    end
    if strcmpi('centered',P1)
        for i = 1:len,
            t_temp = t-(i-1)/Fc;
            ind =  (t_temp >= -1/2/Fc) & (t_temp < 1/2/Fc) ;
            for j1 = 1:size(y,2)   % for each column ...
                x(i,j1) = sum(y(ind,j1))*Fc/Fs;
            end
        end
        x(1,:) = x(1,:)*2;
    elseif strcmpi('left',P1)
        for i = 1:len,
            t_temp = t-(i-1)/Fc;
            ind =  (t_temp >= 0) & (t_temp < 1/Fc) ;
            for j2 = 1:size(y,2)   % for each column ...
                x(i,j2) = sum(y(ind,j2))*Fc/Fs;
            end
        end
    else
        error(message('signal:demod:SignalErr'))
    end
%    w=diff([1; find(diff(y))]);   % <-- a MUCH faster way, but not robust
%    x=w(1:2:length(w))/Fs;
elseif strcmpi(method,'ptm') || strcmpi(method,'ppm'),
    % precondition input by thresholding:
    y = y>.5;
    t = (0:1/Fs:((len-1)/Fs))'*Fc;
    len = ceil( len * Fc / Fs);   % length of message signal
    x = zeros(len,size(y,2));
    for i = 1:len
        t_temp = t-(i-1);
        ind = find( (t_temp >= 0) & (t_temp<1) );
        for j3 = 1:size(y,2)    % for each column ...
            ind2 = y(ind,j3)==1;
            x(i,j3) = t_temp(min(ind(ind2)));
        end
    end
elseif strcmpi(method,'qam'),
    t = (0:1/Fs:((len-1)/Fs))';
    t = t(:,ones(1,size(y,2)));
    x = 2*y.*cos(2*pi*Fc*t);
    x2 = 2*y.*sin(2*pi*Fc*t);
    [b,a]=butter(5,Fc*2/Fs);
    for i = 1:size(y,2),
        x(:,i) = filtfilt(b,a,x(:,i));
        x2(:,i) = filtfilt(b,a,x2(:,i));
    end
    if (r==1),   % convert x2 from a column to a row if necessary
        x2 = x2.';
    end
    
    if isInputSingle
      x2 = single(x2);
    end
end
if (r==1),   % convert x from a column to a row
    x = x.';
end

if isInputSingle
  x = single(x);
end
