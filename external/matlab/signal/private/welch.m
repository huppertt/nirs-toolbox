function varargout = welch(x,esttype,varargin)
%WELCH Welch spectral estimation method.
%   [Pxx,F] = WELCH(X,WINDOW,NOVERLAP,NFFT,Fs,SPECTRUMTYPE,TRACE,ESTTYPE)
%   [Pxx,F] = WELCH({X},WINDOW,NOVERLAP,NFFT,Fs,SPECTRUMTYPE,TRACE,'psd')
%   [Pxx,F] = WELCH({X},WINDOW,NOVERLAP,NFFT,Fs,SPECTRUMTYPE,TRACE,'ms')
%   [Pxy,F] = WELCH({X,Y},WINDOW,NOVERLAP,NFFT,Fs,SPECTRUMTYPE,'cpsd')
%   [Txy,F] = WELCH({X,Y},WINDOW,NOVERLAP,NFFT,Fs,SPECTRUMTYPE,'tfe')
%   [Cxy,F] = WELCH({X,Y},WINDOW,NOVERLAP,NFFT,Fs,SPECTRUMTYPE,'mscohere')
%   [Pxx,F,Pxxc] = WELCH(...)
%   [Pxx,F,Pxxc] = WELCH(...,'ConfidenceLevel',P)
%
%   Inputs:
%      see "help pwelch" for complete description of all input arguments.
%      ESTTYPE - is a string specifying the type of estimate to return, the
%                choices are: psd, cpsd, tfe, and mscohere.
%
%   Outputs:
%      Depends on the input string ESTTYPE:
%      Pxx - Power Spectral Density (PSD) estimate, or
%      MS  - Mean-square spectrum, or
%      Pxy - Cross Power Spectral Density (CPSD) estimate, or
%      Txy - Transfer Function Estimate (TFE), or
%      Cxy - Magnitude Squared Coherence.
%      F   - frequency vector, in Hz if Fs is specified, otherwise it has
%            units of rad/sample

%   Copyright 1988-2014 The MathWorks, Inc.

%   References:
%     [1] Petre Stoica and Randolph Moses, Introduction To Spectral
%         Analysis, Prentice-Hall, 1997, pg. 15
%     [2] Monson Hayes, Statistical Digital Signal Processing and
%         Modeling, John Wiley & Sons, 1996.

narginchk(2,10);
nargoutchk(0,3);

% Parse input arguments.
[x,~,~,y,~,win,winName,winParam,noverlap,k,L,options] = ...
    welchparse(x,esttype,varargin{:});
% Cast to enforce precision rules
options.nfft = signal.internal.sigcasttofloat(options.nfft,'double',...
  'WELCH','NFFT','allownumeric');
noverlap = signal.internal.sigcasttofloat(noverlap,'double','WELCH',...
  'NOVERLAP','allownumeric');
options.Fs = signal.internal.sigcasttofloat(options.Fs,'double','WELCH',...
  'Fs','allownumeric');
k = double(k);

if any([signal.internal.sigcheckfloattype(x,'single')...
    signal.internal.sigcheckfloattype(y,'single'),...
    isa(win,'single')])
  x = single(x);
  y = single(y);
  win = single(win);
end

% Frequency vector was specified, return and plot two-sided PSD
freqVectorSpecified = false; nrow = 1;
if length(options.nfft) > 1,
    freqVectorSpecified = true;
    [~,nrow] = size(options.nfft);
end

% Compute the periodogram power spectrum of each segment and average always
% compute the whole power spectrum, we force Fs = 1 to get a PS not a PSD.

% Initialize
if freqVectorSpecified,
    nFreqs = length(options.nfft);
    if strcmpi(options.range,'onesided')
        warning(message('signal:welch:InconsistentRangeOption'));
    end
    options.range = 'twosided';
else
    nFreqs = options.nfft;
end
% Cast to enforce precision rules
Sxx = zeros(nFreqs,size(x,2),class(x)); %#ok<*ZEROLIKE>

LminusOverlap = L-noverlap;
xStart = 1:LminusOverlap:k*LminusOverlap;
xEnd   = xStart+L-1;

switch esttype,
    case {'ms','power','psd'}
        % Combining method
        if options.maxhold,
            Sxx(:) = -Inf;
            cmethod = @(seg,nextseg) max(seg,k*nextseg); % k will be removed below
        elseif options.minhold,
            Sxx(:) = Inf;
            cmethod = @(seg,nextseg) min(seg,k*nextseg); % k will be removed below
        else
            cmethod = @plus;
        end
        [Pxx,w,units] = localComputeSpectra(Sxx,x,[],xStart,xEnd,win,...
            options,esttype,k,cmethod,freqVectorSpecified);
        
    case 'cpsd'
        numChan = max(size(x,2),size(y,2));
        Sxx = zeros(nFreqs,numChan,class(x)); % Initialize
        cmethod = @plus;
        [Pxx,w,units] = localComputeSpectra(Sxx,x,y,xStart,xEnd,win,...
            options,esttype,k,cmethod,freqVectorSpecified);
      
    case 'tfe'
        numChan = max(size(x,2),size(y,2));
        Sxy = zeros(nFreqs,numChan,class(x));
        cmethod = @plus;
        [Pxx,w,units] = localComputeSpectra(Sxx,x,[],xStart,xEnd,win,...
            options,esttype,k,cmethod,freqVectorSpecified);
        % Cross PSD.  The frequency vector and xunits are not used.
        Pxy = localComputeSpectra(Sxy,y,x,xStart,xEnd,win,...
            options,esttype,k,cmethod,freqVectorSpecified);

        Pxx = bsxfun(@rdivide, Pxy, Pxx);
        
    case 'mscohere'
        % Note: (Sxy1+Sxy2)/(Sxx1+Sxx2) != (Sxy1/Sxy2) + (Sxx1/Sxx2)
        % ie, we can't push the computation of Cxy into computeperiodogram.
        numChan = max(size(x,2),size(y,2));
        Sxy = zeros(nFreqs,numChan,class(x));
        Syy = zeros(nFreqs,size(y,2),class(x));
        cmethod = @plus;
        [Pxx,w,units] = localComputeSpectra(Sxx,x,[],xStart,xEnd,win,...
            options,esttype,k,cmethod,freqVectorSpecified);
        Pyy = localComputeSpectra(Syy,y,[],xStart,xEnd,win,...
            options,esttype,k,cmethod,freqVectorSpecified);
        % Cross PSD.  The frequency vector and xunits are not used.
        Pxy = localComputeSpectra(Sxy,x,y,xStart,xEnd,win,...
            options,esttype,k,cmethod,freqVectorSpecified);

        Pxx = (abs(Pxy).^2)./bsxfun(@times,Pxx,Pyy); % Cxy
end

% Compute confidence intervals if needed.
if ~strcmp(options.conflevel,'omitted')
    if any(strcmpi(esttype,{'ms','power','psd'})) && ~isequal(cmethod,@plus),
        % Always compute the confidence interval around the average spectrum
        cmethod = @plus;
        Sxxc = zeros(nFreqs,size(x,2),class(x));
        [avgPxx,w] = localComputeSpectra(Sxxc,x,[],xStart,xEnd,win,...
            options,esttype,k,cmethod,freqVectorSpecified);
    else
        avgPxx = Pxx;
    end
    Pxxc = confInterval(options.conflevel, avgPxx, x, w, options.Fs, k);
elseif nargout>2
    if any(strcmpi(esttype,{'ms','power','psd'})) && ~isequal(cmethod,@plus),
        % Always compute the confidence interval around the average spectrum
        cmethod = @plus;
        Sxxc = zeros(nFreqs,size(x,2),class(x));
        [avgPxx,w] = localComputeSpectra(Sxxc,x,[],xStart,xEnd,win,...
            options,esttype,k,cmethod,freqVectorSpecified);
    else
        avgPxx = Pxx;
    end
    Pxxc = confInterval(0.95, avgPxx, x, w, options.Fs, k);
else
    Pxxc = [];
end

if nargout==0
    w = {w};
    if strcmpi(units,'Hz'), w = [w,{'Fs',options.Fs}];  end
    % Create a spectrum object to store in the Data object's metadata.
    percOverlap = (noverlap/L)*100;
    hspec = spectrum.welch({winName,winParam},L,percOverlap); %#ok<DWELCH>
    
    switch lower(esttype)
        case 'tfe'
            if strcmpi(options.range,'onesided'), range='half'; else range='whole'; end
            h = dspdata.freqz(Pxx,w{:},'SpectrumRange',range);
        case 'mscohere'
            if strcmpi(options.range,'onesided'), range='half'; else range='whole'; end
            h = dspdata.magnitude(Pxx,w{:},'SpectrumRange',range);
        case 'cpsd'
            h = dspdata.cpsd(Pxx,w{:},'SpectrumType',options.range);
        case {'ms','power'}
            h = dspdata.msspectrum(Pxx,w{:},'SpectrumType',options.range); %#ok<DMSSPEC>
        otherwise
            h = dspdata.psd(Pxx,w{:},'SpectrumType',options.range); %#ok<DPPSD>
    end
    h.Metadata.setsourcespectrum(hspec);
    
    % plot the confidence levels if conflevel is specified.
    if ~isempty(Pxxc)
        h.ConfLevel = options.conflevel;
        h.ConfInterval = Pxxc;
    end
    % center dc component if specified
    if options.centerdc
        centerdc(h);
    end
    plot(h);
    if strcmp(esttype,'power')
        title(getString(message('signal:welch:WelchPowerSpectrumEstimate')));
    end
else
    if options.centerdc
        [Pxx, w, Pxxc] = psdcenterdc(Pxx, w, Pxxc, options);
    end
    
    % If the input is a vector and a row frequency vector was specified,
    % return output as a row vector for backwards compatibility.
    if nrow > 1 && isvector(x)
        Pxx = Pxx.'; w = w.';
    end
    
    % Cast to enforce precision rules   
    % Only cast if output is requested, otherwise, plot using double
    % precision frequency vector.
    if isa(Pxx,'single')
      w = single(w);
    end
    
    if isempty(Pxxc)
        varargout = {Pxx,w}; % Pxx=PSD, MEANSQUARE, CPSD, or TFE
    else
        varargout = {Pxx,w,Pxxc};
    end       
end

end

function [Pxx,w,units] = localComputeSpectra(Sxx,x,y,xStart,xEnd,win,options,esttype,k,cmethod,freqVectorSpecified)

if isempty(y),
    for ii = 1:k
        [Sxxk,w] = computeperiodogram(x(xStart(ii):xEnd(ii),:),win,...
            options.nfft,esttype,options.Fs);
        Sxx  = cmethod(Sxx,Sxxk);
    end
else
    for i = 1:k
        [Sxxk,w] =  computeperiodogram({x(xStart(i):xEnd(i),:),...
            y(xStart(i):xEnd(i),:)},win,options.nfft,esttype,options.Fs);
        Sxx  = cmethod(Sxx,Sxxk);
    end
end
Sxx = Sxx./k; % Average the sum of the periodograms

% Generate the freq vector directly in Hz to avoid roundoff errors due to
% conversions later.
if ~freqVectorSpecified
    w = psdfreqvec('npts',options.nfft, 'Fs',options.Fs);
end

% Compute the 1-sided or 2-sided PSD [Power/freq] or mean-square [Power].
% Also, corresponding freq vector and freq units.
[Pxx,w,units] = computepsd(Sxx,w,options.range,options.nfft,options.Fs,esttype);
end

function Pxxc = confInterval(CL, Pxx, x, w, fs, k)
%   Reference: D.G. Manolakis, V.K. Ingle and S.M. Kagon,
%   Statistical and Adaptive Signal Processing,
%   McGraw-Hill, 2000, Chapter 5
k = fix(k);
c = privatechi2conf(CL,k);

% Cast to enforce precision rules
Pxx = double(Pxx);

PxxcLower = Pxx*c(1);
PxxcUpper = Pxx*c(2);
Pxxc = reshape([PxxcLower; PxxcUpper],size(Pxx,1),2*size(Pxx,2));

% DC and Nyquist bins have only one degree of freedom for real signals
if isreal(x)
    realConf = chi2conf(CL,k/2);
    Pxxc(w == 0,1:2:end) = Pxx(w == 0,:) * realConf(1);
    Pxxc(w == 0,2:2:end) = Pxx(w == 0,:) * realConf(2);
    if isempty(fs)
        Pxxc(w == pi,1:2:end) = Pxx(w == pi,:) * realConf(1);
        Pxxc(w == pi,2:2:end) = Pxx(w == pi,:) * realConf(2);
    else
        Pxxc(w == fs/2,1:2:end) = Pxx(w == fs/2,:) * realConf(1);
        Pxxc(w == fs/2,2:2:end) = Pxx(w == fs/2,:) * realConf(2);
    end
end
end
% [EOF]
