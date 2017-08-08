function [b,a] = fir2(nn, ff, aa, npt, lap, wind)
%FIR2   FIR arbitrary shape filter design using the frequency sampling method.
%   B = FIR2(N,F,A) designs an Nth order linear phase FIR digital filter
%   with the frequency response specified by vectors F and A and returns
%   the filter coefficients in length N+1 vector B.  
%
%   The vectors F and A specify the frequency and magnitude breakpoints for
%   the desired frequency response. The frequencies in F must be given in
%   increasing order with 0.0 < F < 1.0 and 1.0 corresponding to half the
%   sample rate. The first and last elements of F must equal 0 and 1
%   respectively.
%
%   B = FIR2(N,F,A,NPT) specifies the number of points, NPT, for the grid
%   onto which FIR2 linearly interpolates the frequency response. NPT must
%   be greater than 1/2 the filter order (NPT > N/2). If desired, you can
%   interpolate F and A before passing them to FIR2.
%
%   B = FIR2(N,F,A,NPT,LAP) specifies the size of the region, LAP, that
%   FIR2 inserts around duplicate frequency points to provide a smooth but
%   steep transition in the requested frequency response.
%
%   By default FIR2 windows the impulse response with a Hamming window.
%   Other available windows, including Boxcar, Hann, Bartlett, Blackman,
%   Kaiser and Chebwin can be specified with an optional trailing argument.
%
%   For example,
%   B = FIR2(N,F,A,bartlett(N+1)) uses a Bartlett window.
%   B = FIR2(N,F,A,chebwin(N+1,R)) uses a Chebyshev window.
%
%   For filters with a gain other than zero at Fs/2, e.g., highpass
%   and bandstop filters, N must be even.  Otherwise, N will be
%   incremented by one.  In this case the window length should be
%   specified as N+2.
%
%   % Example 1:
%   %   Design a 30th-order lowpass filter and overplot the desired 
%   %   frequency response with the actual frequency response.
%
%   f = [0 0.6 0.6 1];      % Frequency breakpoints
%   m = [1 1 0 0];          % Magnitude breakpoints
%   b = fir2(30,f,m);       % Frequency sampling-based FIR filter design
%   [h,w] = freqz(b,1,128); % Frequency response of filter
%   plot(f,m,w/pi,abs(h))
%   legend('Ideal','fir2 Designed')
%   title('Comparison of Frequency Response Magnitudes')
%
%   % Example 2:
%   %   The chirp.mat file contains a signal, y, that has most of its power 
%   %   above fs/4, or half the Nyquist frequency. Design a 34th-order FIR 
%   %   highpass filter to attenuate the components of the signal below 
%   %   fs/4. Use a cutoff frequency of 0.48.
%
%   load chirp;                     % Load data (y and Fs) into workspace
%   y = y + 0.5*rand(size(y));      % Adding noise
%   f = [0 0.48 0.48 1];            % Frequency breakpoints
%   m = [0 0 1 1];                  % Magnitude breakpoints
%   b = fir2(34,f,m);               % FIR filter design
%   freqz(b,1,512);                 % Frequency response of filter
%   output = filtfilt(b,1,y);       % Zero-phase digital filtering
%   figure;                       
%   subplot(211); plot(y,'b'); title('Original Signal')
%   subplot(212); plot(output,'g'); title('Filtered Signal') 
%
%   See also FIR1, FIRLS, CFIRPM, FIRPM, BUTTER, CHEBY1, CHEBY2, YULEWALK,
%   FREQZ, FILTER, DESIGNFILT.

%   Author(s): L. Shure, 3-27-87
%         C. Denham, 7-26-90, revised
%   Copyright 1988-2013 The MathWorks, Inc.
%     

%   Optional input arguments (see the user's guide):
%     npt - number for interpolation
%     lap - the number of points for jumps in amplitude


narginchk(3,6);

% Cast to enforce precision rules
nn = signal.internal.sigcasttofloat(nn,'double','fir2','N',...
  'allownumeric');

% Check if A is a valid input. 'F' can not be character as the first and
% last elements of F vector should be 0 and 1 respectively so we do not
% need to check its type here. 
aa = signal.internal.sigcasttofloat(aa,'double','fir2','A', 'allownumeric');

[nn,msg1,msg2,msgobj] = firchk(nn,ff(end),aa);

if ~isempty(msg1), error(msgobj); end;
if ~isempty(msg2), warning(msgobj); end;

% Work with filter length instead of filter order
nn = nn + 1;

if (nargin > 3)
  % Cast to enforce precision rules
  npt = signal.internal.sigcasttofloat(npt,'double','fir2','NPT',...
    'allownumeric');
   
    if nargin == 4
        if length(npt) == 1
            if (2 ^ round(log(npt)/log(2))) ~= npt
                % NPT is not an even power of two
                npt = 2^ceil(log(npt)/log(2));
            end
            wind = hamming(nn);
        else
            wind = npt;
            if nn < 1024
                npt = 512;
            else
                npt = 2.^ceil(log(nn)/log(2));
            end
        end
        lap = fix(npt/25);
    elseif nargin == 5
      % Cast to enforce precision rules
      lap = signal.internal.sigcasttofloat(lap,'double','fir2','LAP',...
        'allownumeric');
      
        if length(npt) == 1
            if (2 ^ round(log(npt)/log(2))) ~= npt
                % NPT is not an even power of two
                npt = 2.^ceil(log(npt)/log(2));
            end
            if length(lap) == 1
                wind = hamming(nn);
            else
                wind = lap;
                lap = fix(npt/25);
            end
        else
            wind = npt;
            npt = lap;
            lap = fix(npt/25);
        end
    end
elseif nargin == 3
    if nn < 1024
        npt = 512;
    else
        npt = 2.^ceil(log(nn)/log(2));
    end
    wind = hamming(nn);
    lap = fix(npt/25);
end
% Cast to enforce precision rules
wind = signal.internal.sigcasttofloat(wind,'double','fir2','WINDOW',...
 'allownumeric');

if nn ~= length(wind)
    error(message('signal:fir2:UnequalLengths'))
end

[mf,nf] = size(ff);
[ma,na] = size(aa);
if ma ~= mf || na ~= nf
    error(message('signal:fir2:MismatchedDimensions'))
end
nbrk = max(mf,nf);
if mf < nf
    ff = ff';
    aa = aa';
end

if abs(ff(1)) > eps || abs(ff(nbrk) - 1) > eps
    error(message('signal:fir2:InvalidRange'))
end

ff(1) = 0;
ff(nbrk) = 1;

% interpolate breakpoints onto large grid

H = zeros(1,npt);
nint=nbrk-1;
df = diff(ff'); 
if (any(df < 0))
    error(message('signal:fir2:InvalidFreqVec'))
end

if nn>2*npt
   error(message('signal:fir2:InvalidNpt', ceil( nn/2 ), nn - 1));
end

npt = npt + 1;   % Length of [dc 1 2 ... nyquist] frequencies.

nb = 1;
H(1)=aa(1);
for i=1:nint
    if df(i) == 0
        nb = ceil(nb - lap/2);
        ne = nb + lap;
    else
        ne = fix(ff(i+1)*npt);
    end
    if (nb < 0 || ne > npt)
        error(message('signal:fir2:SignalErr'))
    end
    j=nb:ne;
    if nb == ne
        inc = 0;
    else
        inc = (j-nb)/(ne-nb);
    end
    H(nb:ne) = inc*aa(i+1) + (1 - inc)*aa(i);
    nb = ne + 1;
end

% Fourier time-shift.

dt = 0.5 .* (nn - 1);
rad = -dt .* sqrt(-1) .* pi .* (0:npt-1) ./ (npt-1);
H = H .* exp(rad);

H = [H conj(H(npt-1:-1:2))];   % Fourier transform of real series.
ht = real(ifft(H));            % Symmetric real series.

b = ht(1:nn);         % Raw numerator.
b = b .* wind(:).';   % Apply window.

a = 1;                % Denominator.
