function varargout = spectrogram(x,varargin)
%SPECTROGRAM Spectrogram using a Short-Time Fourier Transform (STFT).
%   S = SPECTROGRAM(X) returns the spectrogram of the signal specified by
%   vector X in the matrix S. By default, X is divided into eight segments
%   with 50% overlap, each segment is windowed with a Hamming window. The
%   number of frequency points used to calculate the discrete Fourier
%   transforms is equal to the maximum of 256 or the next power of two
%   greater than the length of each segment of X.
%
%   If X cannot be divided exactly into eight segments, X will be truncated
%   accordingly.
%
%   S = SPECTROGRAM(X,WINDOW) when WINDOW is a vector, divides X into
%   segments of length equal to the length of WINDOW, and then windows each
%   segment with the vector specified in WINDOW.  If WINDOW is an integer,
%   X is divided into segments of length equal to that integer value, and a
%   Hamming window of equal length is used.  If WINDOW is not specified, the
%   default is used.
%
%   S = SPECTROGRAM(X,WINDOW,NOVERLAP) NOVERLAP is the number of samples
%   each segment of X overlaps. NOVERLAP must be an integer smaller than
%   WINDOW if WINDOW is an integer.  NOVERLAP must be an integer smaller
%   than the length of WINDOW if WINDOW is a vector.  If NOVERLAP is not
%   specified, the default value is used to obtain a 50% overlap.
%
%   S = SPECTROGRAM(X,WINDOW,NOVERLAP,NFFT) specifies the number of
%   frequency points used to calculate the discrete Fourier transforms.
%   If NFFT is not specified, the default NFFT is used. 
%
%   S = SPECTROGRAM(X,WINDOW,NOVERLAP,NFFT,Fs) Fs is the sampling frequency
%   specified in Hz. If Fs is specified as empty, it defaults to 1 Hz. If 
%   it is not specified, normalized frequency is used.
%
%   Each column of S contains an estimate of the short-term, time-localized
%   frequency content of the signal X.  Time increases across the columns
%   of S, from left to right.  Frequency increases down the rows, starting
%   at 0.  If X is a length NX complex signal, S is a complex matrix with
%   NFFT rows and k = fix((NX-NOVERLAP)/(length(WINDOW)-NOVERLAP)) columns.
%   For real X, S has (NFFT/2+1) rows if NFFT is even, and (NFFT+1)/2 rows
%   if NFFT is odd.  
%
%   [S,F,T] = SPECTROGRAM(...) returns a vector of frequencies F and a
%   vector of times T at which the spectrogram is computed. F has length
%   equal to the number of rows of S. T has length k (defined above) and
%   its value corresponds to the center of each segment.
%
%   [S,F,T] = SPECTROGRAM(X,WINDOW,NOVERLAP,F,Fs) computes the two-sided
%   spectrogram at the frequencies specified in vector F.  F must be
%   expressed in hertz and have at least two elements.
%
%   [S,F,T,P] = SPECTROGRAM(...) P is a matrix representing the Power
%   Spectral Density (PSD) of each segment. For real signals, SPECTROGRAM
%   returns the one-sided modified periodogram estimate of the PSD of each
%   segment; for complex signals and in the case when a vector of
%   frequencies is specified, it returns the two-sided PSD.  
%
%   [...]  = SPECTROGRAM(...,SPECTRUMTYPE) uses the window scaling
%   algorithm specified by SPECTRUMTYPE when computing the power spectral
%   density matrix P.
%   SPECTRUMTYPE can be set to 'psd' or 'power':
%     'psd'   - returns the power spectral density
%     'power' - scales each estimate of the PSD by the equivalent noise
%               bandwidth of the window (in hertz).  Use this option to
%               obtain an estimate of the power at each frequency.
%   The default value for SPECTRUMTYPE is 'psd'.
%
%   [...] = SPECTROGRAM(...,FREQRANGE)  returns the PSD over the specified
%   range of frequencies based upon the value of FREQRANGE:
%
%      'onesided' - returns the one-sided matrix P of a real input signal X.
%         If NFFT is even, P has length NFFT/2+1 and is computed over the
%         interval [0,pi].  If NFFT is odd, P has length (NFFT+1)/2 and
%         is computed over the interval [0,pi). When Fs is specified, the
%         intervals become [0,Fs/2) and [0,Fs/2] for even and odd NFFT,
%         respectively.
%
%      'twosided' - returns the two-sided matrix P for either real or complex
%         input X.  P has length NFFT and is computed over the interval
%         [0,2*pi). When Fs is specified, the interval becomes [0,Fs).
%
%      'centered' - returns the centered two-sided matrix P for either real 
%         or complex X.  P has length NFFT and is computed over the interval
%         (-pi, pi] for even length NFFT and (-pi, pi) for odd length NFFT.
%         When Fs is specified, the intervals become (-Fs/2, Fs/2] and
%         (-Fs/2, Fs/2) for even and odd NFFT, respectively.
%
%      FREQRANGE may be placed in any position in the input argument list
%      after NOVERLAP.  The default value of FREQRANGE is 'onesided' when X
%      is real and 'twosided' when X is complex.
%
%   SPECTROGRAM(...) with no output arguments plots the PSD estimate for
%   each segment on a surface in the current figure. It uses
%   SURF(f,t,10*log10(abs(P)) where P is the fourth output argument. A
%   trailing input string, FREQLOCATION, controls where MATLAB displays the
%   frequency axis. This string can be either 'xaxis' or 'yaxis'.  Setting
%   this FREQLOCATION to 'yaxis' displays frequency on the y-axis and time
%   on the x-axis.  The default is 'xaxis' which displays the frequency on
%   the x-axis. If FREQLOCATION is specified when output arguments are
%   requested, it is ignored.
%
%   EXAMPLE 1: Display the PSD of each segment of a quadratic chirp.
%     t=0:0.001:2;                    % 2 secs @ 1kHz sample rate
%     y=chirp(t,100,1,200,'q');       % Start @ 100Hz, cross 200Hz at t=1sec 
%     spectrogram(y,128,120,128,1E3); % Display the spectrogram
%     title('Quadratic Chirp: start at 100Hz and cross 200Hz at t=1sec');
%
%   EXAMPLE 2: Waterfall display of the PSD of each segment of a VCO
%     Fs = 10e3;
%     t = 0:1/Fs:2;
%     x1 = vco(sawtooth(2*pi*t,0.5),[0.1 0.4]*Fs,Fs);
%     spectrogram(x1,kaiser(256,5),220,512,Fs,'yaxis');
%     view(-45,65)
%     colormap bone
%
%   See also PERIODOGRAM, PWELCH, SPECTRUM, GOERTZEL.

% [1] Oppenheim, A.V., and R.W. Schafer, Discrete-Time Signal Processing,
% Prentice-Hall, Englewood Cliffs, NJ, 1989, pp. 713-718.
% [2] Mitra, S. K., Digital Signal Processing. A Computer-Based Approach.
% 2nd Ed. McGraw-Hill, N.Y., 2001.

%   Copyright 1988-2014 The MathWorks, Inc.

narginchk(1,8);
nargoutchk(0,4);

% Frequency axis location flag; Handled up-front so that we can remove it
% from varargin and then reuse welchparse. Ignored in case when outputs are
% requested.
faxisloc = 'xaxis';
if (nargin > 1 && ischar(varargin{end})) && any(strcmpi(varargin{end},{'yaxis','xaxis'})),
    if strcmpi(varargin{end},'yaxis'),
        faxisloc = 'yaxis';
    end
    varargin(end)=[];
end

% Check for valid input signal
chkinput(x);

% Parse input arguments (using the PWELCH parser since we share the same API).
% The following outputs are NOT used: y, Ly, winName,winParam, k, and L. 
[esttype, varargin] = psdesttype({'psd','power'},'psd',varargin); % Look for psd and power flags
[x,nx,xisreal,y,Ly,win,winName,winParam,noverlap,k,L,options] = ...
    welchparse(x,esttype,varargin{:});  %#ok

% Determine whether an empty was specified for Fs (i.e., Fs=1Hz) or
% returned by welchparse which means normalized Fs is used.
Fs = options.Fs; isFsnormalized = false;

% Cast to enforce Precision rules
Fs = double(Fs);
noverlap = signal.internal.sigcasttofloat(noverlap,'double',...
  'spectrogram','NOVERLAP','allownumeric');

if isempty(Fs)
    lv = length(varargin);
    if lv == 5,
        isFsnormalized = false; % Fs = 1Hz.
        Fs = 1;
    elseif lv < 5
        isFsnormalized = true;  % Fs is normalized
        Fs = 2*pi;
    end
end

% Window length
nwind = length(win);

% Make x and win into columns
x = x(:); 
win = win(:); 

% Determine the number of columns of the STFT output (i.e., the S output)
ncol = fix((nx-noverlap)/(nwind-noverlap));

%
% Pre-process X
%
colindex = 1 + (0:(ncol-1))*(nwind-noverlap);
rowindex = (1:nwind)';
% 'xin' should be of the same datatype as 'x'
xin = zeros(nwind,ncol,class(x)); %#ok<*ZEROLIKE>

% Put x into columns of xin with the proper offset
xin(:) = x(rowindex(:,ones(1,ncol))+colindex(ones(nwind,1),:)-1);

% Apply the window to the array of offset signal segments.
xin = win(:,ones(1,ncol)).*xin;

%
% Compute the raw STFT with the appropriate algorithm
%
% Cast to enforce Precision rules
nfft = double(options.nfft);
 
freqvecspecified = false;
if length(nfft) > 1,
    % Frequency vector was specified, return and plot two-sided PSD
    freqvecspecified = true; 
    if strcmpi(options.range,'onesided')
        warning(message('signal:welch:InconsistentRangeOption'));
    end
    options.range = 'twosided';
end

[y,f] = computeDFT(xin,nfft,Fs);

% Cast to enforce precision rules
if (isa(xin,'single'))
  y = single(y);
end

if ~freqvecspecified && strcmpi(options.range,'onesided')
    f = psdfreqvec('npts',nfft,'Fs',Fs,'Range','half');
    y = y(1:length(f),:);
end

% colindex already takes into account the noverlap factor; Return a T
% vector whose elements are centered in the segment.
t = ((colindex-1)+((nwind)/2)')/Fs; 

% Outputs
switch nargout,
    case 0
        % Use surface to display spectrogram
        [Pxx,f] = compute_PSD(win,y,nfft,f,options,Fs,esttype); 
        
        % Surf encodes data points at the edges and takes the color from
        % the last edge so we need to add an additional point so that surf
        % does the right thing. This is important especially when
        % spectrogram has only one estimate (e.g. window length = signal
        % length). Although this issue also exists on the frequency
        % direction we will not add an extra row since we will never
        % encounter a single frequency point estimate. For the plot we set
        % time values to be at: nwin/2-a/2, nwin/2+a/2, nwin/2+3a/2,
        % nwin/2+5a/2 ... where a is the number of new samples for each
        % segment (i.e. nwind-noverlap). For the case of zero overlap this
        % corresponds to 0, nwind, 2*nwind, ...        
        a = nwind - noverlap;
        t = [(nwind/2-a/2)/Fs  t+((a/2)/Fs)];                       
        Pxx = [Pxx Pxx(:,end)];         
        
        displayspectrogram(t,f,double(Pxx),isFsnormalized,faxisloc,esttype);           
        
    case 1
        if options.centerdc && length(options.nfft)==1
            y = psdcenterdc(y,f,[],options);
        end
        varargout = {y};
    case 2
        if options.centerdc && length(options.nfft)==1
            [y,f] = psdcenterdc(y,f,[],options);
        end
        % Cast to enforce precision rules.
        if isa(y,'single')
          f = single(f);
        end
        varargout = {y,f};
    case 3
        if options.centerdc && length(options.nfft)==1
            [y,f] = psdcenterdc(y,f,[],options);
        end
        % Cast to enforce precision rules.      
        if isa(y,'single')
          f = single(f);
          t = single(t);
        end      
        varargout = {y,f,t};
    case 4      
        [Pxx,f] = compute_PSD(win,y,nfft,f,options,Fs,esttype);
        % Cast to enforce precision rules.        
        if isa(Pxx,'single')
          f = single(f);
          t = single(t);
        end              
        varargout = {y,f,t,Pxx};
end


%--------------------------------------------------------------------------
function chkinput(x)
% Check for valid input signal

if isempty(x) || issparse(x) || ~isfloat(x),
    error(message('signal:spectrogram:MustBeFloat', 'X'));
end

if min(size(x))~=1,
    error(message('signal:spectrogram:MustBeVector', 'X'));
end


%--------------------------------------------------------------------------
function displayspectrogram(t,f,Pxx,isFsnormalized,faxisloc,esttype)

% Cell array of the standard frequency units strings

if isFsnormalized, 
    f = f/pi; % Normalize the freq axis
    frequnitstrs = getfrequnitstrs;
    freqlbl = frequnitstrs{1};
else
    [f,~,uf] = engunits(f,'unicode');
    freqlbl = getfreqlbl([uf 'Hz']);
end

% Use engineering units
[t,~,ut] = engunits(t,'unicode','time');

newplot;
if strcmpi(faxisloc,'yaxis'),    
    args = {t,f,10*log10(abs(Pxx)+eps)};    
    % Axis labels
    xlbl = [getString(message('signal:spectrogram:Time')) ' (' ut ')'];
    ylbl = freqlbl;
else
    args = {f,t,10*log10(abs(Pxx')+eps)};
    xlbl = freqlbl;
    ylbl = [getString(message('signal:spectrogram:Time')) ' (' ut ')'];
end
hndl = surf(args{:},'EdgeColor','none'); %#ok<NASGU>
% AZ = 0, EL = 90 is directly overhead and the default 2-D view.
view(0,90);

axis xy; axis tight;
colormap('default')
if strcmpi(esttype,'power')
    cblabel = getString(message('signal:dspdata:dspdata:PowerdB'));
else
    if isFsnormalized
        cblabel = getString(message('signal:dspdata:dspdata:PowerfrequencydBradsample'));
    else
        cblabel = getString(message('signal:dspdata:dspdata:PowerfrequencydBHz'));
    end
end
%sigutils.internal.colorbari('titlelong',cblabel);
h = colorbar;
h.Label.String = cblabel;

ylabel(ylbl);
xlabel(xlbl);

% -------------------------------------------------------------------------
function [Pxx,W] = compute_PSD(win,y,nfft,f,options,Fs,esttype)
% Evaluate the window normalization constant.  A 1/N factor has been
% omitted since it will cancel below.
if strcmpi(esttype,'power')
    % The window is convolved with every power spectrum peak, therefore
    % compensate for the DC value squared to obtain correct peak heights.
    U = sum(win)^2;
else
    U = win'*win;  % compensates for the power of the window.
end
Sxx = y.*conj(y)/U; % Auto spectrum.

% The computepsd function expects NFFT to be a scalar
if length(nfft) > 1, nfft = length(nfft); end

% Compute the one-sided or two-sided PSD [Power/freq]. Also compute
% the corresponding half or whole power spectrum [Power].
[Pxx,W] = computepsd(Sxx,f,options.range,nfft,Fs,esttype);

if options.centerdc && length(options.nfft)==1
    [Pxx, W] = psdcenterdc(Pxx, W, [], options);
end


% [EOF]
