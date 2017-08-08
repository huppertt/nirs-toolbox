function [Spec,f] = spectrum(varargin)
%SPECTRUM Spectral Estimation.
%
%   SPECTRUM will be removed in a future release of MATLAB.  Use one of the
%   following functions instead:
%      <a href="matlab:help periodogram">periodogram</a>
%      <a href="matlab:help pwelch">pwelch</a>
%      <a href="matlab:help pburg">pburg</a>
%      <a href="matlab:help pcov">pcov</a>
%      <a href="matlab:help pmcov">pmcov</a>
%      <a href="matlab:help pyulear">pyulear</a>
%      <a href="matlab:help pmtm">pmtm</a>
%      <a href="matlab:help peig">peig</a>
%      <a href="matlab:help pmusic">pmusic</a>
%      <a href="matlab:help tfestimate">tfestimate</a>
%      <a href="matlab:help mscohere">mscohere</a>
%
%   <a href="matlab:help signal">Signal Processing Toolbox TOC</a>

%   H = SPECTRUM.<ESTIMATOR> returns a spectral estimator in H of type
%   specified by ESTIMATOR.
%
%   There are three types of estimators: power spectral density (PSD),
%   mean-square spectrum (MSS), and pseudo spectrum.  Valid values for
%   <ESTIMATOR> for each of these three types are listed below.   Type
%   "help spectrum.<ESTIMATOR>" to get help on a specific spectrum
%   estimator e.g., "help spectrum.welch":
%
%   ----------------------------
%   Power spectral density (PSD)
%   ----------------------------
%   <a href="matlab:help spectrum.burg">burg</a> - Burg                             <a href="matlab:help spectrum.periodogram">periodogram</a> - Periodogram
%   <a href="matlab:help spectrum.cov">cov</a>  - Covariance                       <a href="matlab:help spectrum.welch">welch</a>       - Welch 
%   <a href="matlab:help spectrum.mcov">mcov</a> - Modified Covariance              <a href="matlab:help spectrum.yulear">yulear</a>      - Yule-Walker AR
%   <a href="matlab:help spectrum.mtm">mtm</a>  - Thomson multitaper method (MTM)
%
%   To calculate the PSD you must create one of the estimators above and
%   pass it to one of the functions below along with the data.
%
%   <a href="matlab:help spectrum/psd">psd</a>     - calculates the PSD
%   <a href="matlab:help spectrum/psdopts">psdopts</a> - options to calculate the PSD
%
%   --------------------------
%   Mean-square spectrum (MSS)                        
%   --------------------------
%   <a href="matlab:help spectrum.periodogram">periodogram</a> - Periodogram 
%   <a href="matlab:help spectrum.welch">welch</a>       - Welch
%
%   To calculate the MSS you must create one of the estimators above and
%   pass it to one of the functions below along with the data.
%
%   <a href="matlab:help spectrum/msspectrum">msspectrum</a>     - calculates the Mean-squared Spectrum (MSS)
%   <a href="matlab:help spectrum/msspectrumopts">msspectrumopts</a> - options to calculate the MSS
%
%   ---------------
%   Pseudo spectrum
%   ---------------
%   <a href="matlab:help spectrum.eigenvector">eigenvector</a> - Eigenvector
%   <a href="matlab:help spectrum.music">music</a>       - Multiple Signal Classification
%
%   To calculate the pseudo spectrum you must create one of the estimators
%   above and pass it to one of the functions below along with the data.
%
%   <a href="matlab:help spectrum/powerest">powerest</a>           - computes the powers and frequencies of sinusoids
%   <a href="matlab:help spectrum/pseudospectrum">pseudospectrum</a>     - calculates the pseudospectrum
%   <a href="matlab:help spectrum/pseudospectrumopts.">pseudospectrumopts</a> - options to calculate the pseudospectrum
%
%   EXAMPLE: Spectral analysis of a signal that contains a 200Hz cosine
%            % plus noise.
%            Fs = 1000;   t = 0:1/Fs:.296;
%            x = cos(2*pi*t*200)+randn(size(t));  
%            h = spectrum.welch;    % Create a Welch spectral estimator. 
%
%            Hpsd = psd(h,x,'Fs',Fs);             % Calculate the PSD 
%            plot(Hpsd)                           % Plot the PSD.
%
%
%   See also DSPDATA.

%   Author(s): J.N. Little, 7-9-86
%   	   C. Denham, 4-25-88, revised
%   	   L. Shure, 12-20-88, revised
%   	   J.N. Little, 8-31-89, revised 
%   	   L. Shure, 8-11-92, revised 
%   	   T. Krauss, 4-15-93, revised
%   Copyright 1988-2012 The MathWorks, Inc.

%   The units on the power spectra Pxx and Pyy are such that, using
%   Parseval's theorem: 
%
%        SUM(Pxx)/LENGTH(Pxx) = SUM(X.^2)/LENGTH(X) = COV(X)
%
%   The RMS value of the signal is the square root of this.
%   If the input signal is in Volts as a function of time, then
%   the units on Pxx are Volts^2*seconds = Volt^2/Hz.
%
%   Here are the covariance, RMS, and spectral amplitude values of
%   some common functions:
%         Function   Cov=SUM(Pxx)/LENGTH(Pxx)   RMS        Pxx
%         a*sin(w*t)        a^2/2            a/sqrt(2)   a^2*LENGTH(Pxx)/4
%Normal:  a*rand(t)         a^2              a           a^2
%Uniform: a*rand(t)         a^2/12           a/sqrt(12)  a^2/12
%   
%   For example, a pure sine wave with amplitude A has an RMS value
%   of A/sqrt(2), so A = SQRT(2*SUM(Pxx)/LENGTH(Pxx)).
%
%   See Page 556, A.V. Oppenheim and R.W. Schafer, Digital Signal
%   Processing, Prentice-Hall, 1975.

% Obsolete in 2007a.
warning(message('signal:spectrum:obsoleteFunction'));

narginchk(1,8)
[x,y,nfft,noverlap,window,Fs,p,dflag]=specchk(varargin);

if isempty(p),
	p = .95;   % default confidence interval even if not asked for
end

n = length(x);		% Number of data points
nwind = length(window);
if n < nwind    % zero-pad x (and y) if length less than the window length
    x(nwind)=0;  n=nwind;
    if ~isempty(y), y(nwind)=0;  end
end
x = x(:);		% Make sure x and y are column vectors
y = y(:);
k = fix((n-noverlap)/(nwind-noverlap));	% Number of windows
					% (k = fix(n/nwind) for noverlap=0)
index = 1:nwind;
KMU = k*norm(window)^2;	% Normalizing scale factor ==> asymptotically unbiased
% KMU = k*sum(window)^2;% alt. Nrmlzng scale factor ==> peaks are about right

if (isempty(y))	% Single sequence case.
	Pxx = zeros(nfft,1); Pxx2 = zeros(nfft,1);
	for i=1:k
                if strcmp(dflag,'linear')
                    xw = window.*detrend(x(index));
                elseif strcmp(dflag,'none')
                    xw = window.*(x(index));
                else
                    xw = window.*detrend(x(index),0);
                end
		index = index + (nwind - noverlap);
		Xx = abs(fft(xw,nfft)).^2;
		Pxx = Pxx + Xx;
		Pxx2 = Pxx2 + abs(Xx).^2;
	end
	% Select first half
	if ~any(any(imag(x)~=0)),   % if x and y are not complex
		if rem(nfft,2),    % nfft odd
			select = 1:(nfft+1)/2;
		else
			select = 1:nfft/2+1;   % include DC AND Nyquist
		end
	else
		select = 1:nfft;
	end
	Pxx = Pxx(select);
	Pxx2 = Pxx2(select);
	cPxx = zeros(size(Pxx));
	if k > 1
		c = (k.*Pxx2-abs(Pxx).^2)./(k-1);
		c = max(c,zeros(size(Pxx)));
		cPxx = sqrt(c);
	end
	ff = sqrt(2)*erfinv(p);  % Equal-tails.
	Pxx = Pxx/KMU;
	Pxxc = ff.*cPxx/KMU;
	P = [Pxx Pxxc];
else
	Pxx = zeros(nfft,1); % Dual sequence case.
	Pyy = Pxx; Pxy = Pxx; Pxx2 = Pxx; Pyy2 = Pxx; Pxy2 = Pxx;

	for i=1:k
                if strcmp(dflag,'linear')
                    xw = window.*detrend(x(index));
		    yw = window.*detrend(y(index));
                elseif strcmp(dflag,'none')
                    xw = window.*(x(index));
		    yw = window.*(y(index));
                else
                    xw = window.*detrend(x(index),0);
		    yw = window.*detrend(y(index),0);
                end
		index = index + (nwind - noverlap);
		Xx = fft(xw,nfft);
		Yy = fft(yw,nfft);
		Yy2 = abs(Yy).^2;
		Xx2 = abs(Xx).^2;
		Xy  = Yy .* conj(Xx);
		Pxx = Pxx + Xx2;
		Pyy = Pyy + Yy2;
		Pxy = Pxy + Xy;
		Pxx2 = Pxx2 + abs(Xx2).^2;
		Pyy2 = Pyy2 + abs(Yy2).^2;
		Pxy2 = Pxy2 + Xy .* conj(Xy);
	end

	% Select first half
	if ~any(any(imag([x y])~=0)),   % if x and y are not complex
		if rem(nfft,2),    % nfft odd
			select = 1:(nfft+1)/2;
		else
			select = 1:nfft/2+1;   % include DC AND Nyquist
		end
	else
		select = 1:nfft;
	end

	Pxx = Pxx(select);
	Pyy = Pyy(select);
	Pxy = Pxy(select);
	Pxx2 = Pxx2(select);
	Pyy2 = Pyy2(select);
	Pxy2 = Pxy2(select);
	
	cPxx = zeros(size(Pxx));
	cPyy = cPxx;
	cPxy = cPxx;
	if k > 1
   		c = max((k.*Pxx2-abs(Pxx).^2)./(k-1),zeros(size(Pxx)));
   		cPxx = sqrt(c);
   		c = max((k.*Pyy2-abs(Pyy).^2)./(k-1),zeros(size(Pxx)));
   		cPyy = sqrt(c);
   		c = max((k.*Pxy2-abs(Pxy).^2)./(k-1),zeros(size(Pxx)));
   		cPxy = sqrt(c);
	end

	Txy = Pxy./Pxx;
	Cxy = (abs(Pxy).^2)./(Pxx.*Pyy);
	
	ff = sqrt(2)*erfinv(p);  % Equal-tails.

	Pxx = Pxx/KMU;
	Pyy = Pyy/KMU;
	Pxy = Pxy/KMU;
	Pxxc = ff.*cPxx/KMU;
	Pxyc = ff.*cPxy/KMU;
	Pyyc = ff.*cPyy/KMU;
	P = [Pxx Pyy Pxy Txy Cxy Pxxc Pyyc Pxyc];
end

freq_vector = (select - 1)'*Fs/nfft;
if nargout == 0,   % do plots
        newplot;
	c = [max(Pxx-Pxxc,0)  Pxx+Pxxc];
	c = c.*(c>0);
	semilogy(freq_vector,Pxx,freq_vector,c(:,1),'--',...
		freq_vector,c(:,2),'--');
	title('Pxx - X Power Spectral Density')
	xlabel('Frequency')
	if (isempty(y)),   % single sequence case
		return
	end
	pause

        newplot;
	c = [max(Pyy-Pyyc,0)  Pyy+Pyyc];
	c = c.*(c>0);
	semilogy(freq_vector,Pyy,freq_vector,c(:,1),'--',...
		freq_vector,c(:,2),'--');
	title('Pyy - Y Power Spectral Density')
	xlabel('Frequency')
	pause

        newplot;
	semilogy(freq_vector,abs(Txy));
	title('Txy - Transfer function magnitude')
	xlabel('Frequency')
	pause

        newplot;
	plot(freq_vector,180/pi*angle(Txy)), ...
	title('Txy - Transfer function phase')
	xlabel('Frequency')
	pause

        newplot;
	plot(freq_vector,Cxy);
 	title('Cxy - Coherence')
 	xlabel('Frequency')

elseif nargout ==1, 
	Spec = P;
elseif nargout ==2, 
	Spec = P;
	f = freq_vector;
end

function [x,y,nfft,noverlap,window,Fs,p,dflag] = specchk(P)
%SPECCHK Helper function for SPECTRUM
%   SPECCHK(P) takes the cell array P and uses each cell as 
%   an input argument.  Assumes P has between 1 and 7 elements.

%   Author(s): T. Krauss, 4-6-93

if length(P{1})<=1
    error(message('signal:spectrum:InvalidInput'));
elseif (length(P)>1),
    if (all(size(P{1})==size(P{2})) && (length(P{1})>1) ) || ...
       length(P{2})>1,   % 0ne signal or 2 present?
        % two signals, x and y, present
        x = P{1}; y = P{2}; 
        % shift parameters one left
        P(1) = [];
    else 
        % only one signal, x, present
        x = P{1}; y = []; 
    end
else  % length(P) == 1
    % only one signal, x, present
    x = P{1}; y = []; 
end

% now x and y are defined; let's get the rest

if length(P) == 1
    nfft = min(length(x),256);
    window = hanning(nfft);
    noverlap = 0;
    Fs = 2;
    p = [];
    dflag = 'linear';
elseif length(P) == 2
    if isempty(P{2}),    dflag = 'linear'; nfft = min(length(x),256); 
    elseif ischar(P{2}), dflag = P{2};     nfft = min(length(x),256); 
    else              dflag = 'linear'; nfft = P{2};   end
    window = hanning(nfft);
    noverlap = 0;
    Fs = 2;
    p = [];
elseif length(P) == 3
    if isempty(P{2}), nfft = min(length(x),256); else nfft=P{2};     end
    if isempty(P{3}),    dflag = 'linear'; noverlap = 0;
    elseif ischar(P{3}), dflag = P{3};     noverlap = 0;
    else              dflag = 'linear'; noverlap = P{3}; end
    window = hanning(nfft);
    Fs = 2;
    p = [];
elseif length(P) == 4
    if isempty(P{2}), nfft = min(length(x),256); else nfft=P{2};     end
    if ischar(P{4})
        dflag = P{4};
        window = hanning(nfft);
    else
        dflag = 'linear';
        window = P{4};  window = window(:);   % force window to be a column
        if length(window) == 1, window = hanning(window); end
        if isempty(window), window = hanning(nfft); end
    end
    if isempty(P{3}), noverlap = 0;  else noverlap=P{3}; end
    Fs = 2;
    p = [];
elseif length(P) == 5
    if isempty(P{2}), nfft = min(length(x),256); else nfft=P{2};     end
    window = P{4};  window = window(:);   % force window to be a column
    if length(window) == 1, window = hanning(window); end
    if isempty(window), window = hanning(nfft); end
    if isempty(P{3}), noverlap = 0;  else noverlap=P{3}; end
    if ischar(P{5})
        dflag = P{5};
        Fs = 2;
    else
        dflag = 'linear';
        if isempty(P{5}), Fs = 2; else Fs = P{5}; end
    end
    p = [];
elseif length(P) == 6
    if isempty(P{2}), nfft = min(length(x),256); else nfft=P{2};     end
    window = P{4};  window = window(:);   % force window to be a column
    if length(window) == 1, window = hanning(window); end
    if isempty(window), window = hanning(nfft); end
    if isempty(P{3}), noverlap = 0;  else noverlap=P{3}; end
    if isempty(P{5}), Fs = 2;     else    Fs = P{5}; end
    if ischar(P{6})
        dflag = P{6};
        p = [];
    else
        dflag = 'linear';
        if isempty(P{6}), p = .95;    else    p = P{6}; end
    end
elseif length(P) == 7
    if isempty(P{2}), nfft = min(length(x),256); else nfft=P{2};     end
    window = P{4};  window = window(:);   % force window to be a column
    if length(window) == 1, window = hanning(window); end
    if isempty(window), window = hanning(nfft); end
    if isempty(P{3}), noverlap = 0;  else noverlap=P{3}; end
    if isempty(P{5}), Fs = 2;     else    Fs = P{5}; end
    if isempty(P{6}), p = .95;    else    p = P{6}; end
    if ischar(P{7})
        dflag = P{7};
    else
        error(message('signal:spectrum:NeedStringDFLAG'));
    end
end

% NOW do error checking
if (nfft<length(window)), 
    error(message('signal:spectrum:WindowTooBig'));
end
if (noverlap >= length(window)),
    error(message('signal:spectrum:NOVERLAPTooBig'));
end
if (nfft ~= abs(round(nfft)))||(noverlap ~= abs(round(noverlap))),
    error(message('signal:spectrum:NeedPosInts'));
end
if ~isempty(p),
    if (numel(p)>1)||(p(1,1)>1)||(p(1,1)<0),
        error(message('signal:spectrum:StrictlyInsideUnitRange'));
    end
end
if min(size(x))~=1,
    error(message('signal:spectrum:NeedVecInput'));
end
if (min(size(y))~=1)&&(~isempty(y)),
    error(message('signal:spectrum:NeedVecInput'));
end
if (length(x)~=length(y))&&(~isempty(y)),
    error(message('signal:spectrum:MismatchedDimensions'));
end
