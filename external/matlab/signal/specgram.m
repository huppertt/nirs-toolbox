function [yo,fo,to] = specgram(varargin)
%SPECGRAM Spectrogram using a Short-Time Fourier Transform (STFT).
%   SPECGRAM has been replaced by SPECTROGRAM.  SPECGRAM still works but
%   may be removed in the future. Use SPECTROGRAM instead. Type help
%   SPECTROGRAM for details.
%
%   See also PERIODOGRAM, SPECTRUM/PERIODOGRAM, PWELCH, SPECTRUM/WELCH, GOERTZEL.

%   Author(s): L. Shure, 1-1-91
%              T. Krauss, 4-2-93, updated
%   Copyright 1988-2010 The MathWorks, Inc.

error(nargchk(1,5,nargin,'struct'))
[x,nfft,Fs,window,noverlap]=specgramchk(varargin);
    
nx = length(x);
nwind = length(window);
if nx < nwind    % zero-pad x if it has length less than the window length
    x(nwind)=0;  nx=nwind;
end
x = x(:); % make a column vector for ease later
window = window(:); % be consistent with data set

ncol = fix((nx-noverlap)/(nwind-noverlap));
colindex = 1 + (0:(ncol-1))*(nwind-noverlap);
rowindex = (1:nwind)';
if length(x)<(nwind+colindex(ncol)-1)
    x(nwind+colindex(ncol)-1) = 0;   % zero-pad x
end

if length(nfft)>1
    df = diff(nfft);
    evenly_spaced = all(abs(df-df(1))/Fs<1e-12);  % evenly spaced flag (boolean)
    use_chirp = evenly_spaced & (length(nfft)>20);
else
    use_chirp = 0;
end

if (length(nfft)==1) || use_chirp
    y = zeros(nwind,ncol);

    % put x into columns of y with the proper offset
    % should be able to do this with fancy indexing!
    y(:) = x(rowindex(:,ones(1,ncol))+colindex(ones(nwind,1),:)-1);

    % Apply the window to the array of offset signal segments.
    y = window(:,ones(1,ncol)).*y;

    if ~use_chirp     % USE FFT
        % now fft y which does the columns
        y = fft(y,nfft);
        if ~any(any(imag(x)))    % x purely real
            if rem(nfft,2),    % nfft odd
                select = 1:(nfft+1)/2;
            else
                select = 1:nfft/2+1;
            end
            y = y(select,:);
        else
            select = 1:nfft;
        end
        f = (select - 1)'*Fs/nfft;
    else % USE CHIRP Z TRANSFORM
        f = nfft(:);
        f1 = f(1);
        f2 = f(end);
        m = length(f);
        w = exp(-1i*2*pi*(f2-f1)/(m*Fs));
        a = exp(1i*2*pi*f1/Fs);
        y = czt(y,m,w,a);
    end
else  % evaluate DFT on given set of frequencies
    f = nfft(:);
    q = nwind - noverlap;
    extras = floor(nwind/q);
    x = [zeros(q-rem(nwind,q)+1,1); x];
    % create windowed DTFT matrix (filter bank)
    D = window(:,ones(1,length(f))).*exp((-1i*2*pi/Fs*((nwind-1):-1:0)).'*f'); 
    y = upfirdn(x,D,1,q).';
    y(:,[1:extras+1 end-extras+1:end]) = []; 
end

t = (colindex-1)'/Fs;

% take abs, and use image to display results
if nargout == 0
    newplot;
    if length(t)==1
        imagesc([0 1/f(2)],f,20*log10(abs(y)+eps));axis xy; colormap(jet)
    else
        % Shift time vector by half window length; the overlap factor has
        % already been accounted for in the colindex variable.
        t = ((colindex-1)+((nwind)/2)')/Fs; 
        imagesc(t,f,20*log10(abs(y)+eps));axis xy; colormap(jet)
    end
    xlabel('Time')
    ylabel('Frequency')
elseif nargout == 1,
    yo = y;
elseif nargout == 2,
    yo = y;
    fo = f;
elseif nargout == 3,
    yo = y;
    fo = f;
    to = t;
end

function [x,nfft,Fs,window,noverlap] = specgramchk(P)
%SPECGRAMCHK Helper function for SPECGRAM.
%   SPECGRAMCHK(P) takes the cell array P and uses each cell as 
%   an input argument.  Assumes P has between 1 and 5 elements.

x = P{1}; 
if (length(P) > 1) && ~isempty(P{2})
    nfft = P{2};
else
    nfft = min(length(x),256);
end
if (length(P) > 2) && ~isempty(P{3})
    Fs = P{3};
else
    Fs = 2;
end
if length(P) > 3 && ~isempty(P{4})
    window = P{4}; 
else
    if length(nfft) == 1
        window = hanning(nfft);
    else
      error(message('signal:specgram:NeedWindow'));
    end
end
if length(window) == 1, window = hanning(window); end
if (length(P) > 4) && ~isempty(P{5})
    noverlap = P{5};
else
    noverlap = ceil(length(window)/2);
end

% NOW do error checking
if (length(nfft)==1) && (nfft<length(window)), 
    error(message('signal:specgram:WindowTooBig'));
end
if (noverlap >= length(window)),
    error(message('signal:specgram:OverlapTooBig'));
end
if (length(nfft)==1) && (nfft ~= abs(round(nfft)))
    error(message('signal:specgram:NFFTMustBePositive'));
end
if (noverlap ~= abs(round(noverlap))),
    error(message('signal:specgram:OverlapMustBePositive'));
end
if min(size(x))~=1,
    error(message('signal:specgram:MustBeVector'));
end

