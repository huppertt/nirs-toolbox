function varargout = pmtm(x,varargin)
%PMTM   Power Spectral Density (PSD) estimate via the Thomson multitaper 
%   method (MTM).
%   Pxx = PMTM(X) returns the Power Spectral Density (PSD) estimate, Pxx,
%   of a discrete-time signal X. When X is a vector, it is converted to a
%   column vector and treated as a single channel.  When X is a matrix, the
%   PSD is computed independently for each column and stored in the
%   corresponding column of Pxx. Pxx is the distribution of power per unit
%   frequency. The frequency is expressed in units of radians/sample.
%
%   For real signals, PMTM returns the one-sided PSD by default; for 
%   complex signals, it returns the two-sided PSD.  Note that a one-sided 
%   PSD contains the total power of the input signal.
%
%   Pxx = PMTM(X,NW) specifies NW as the "time-bandwidth product" for the
%   discrete prolate spheroidal sequences (or Slepian sequences) used as 
%   data windows.  Typical choices for NW are 2, 5/2, 3, 7/2, or 4. If
%   empty or omitted, NW defaults to 4. By default, PMTM drops the last
%   taper because its corresponding eigenvalue is significantly smaller
%   than 1. Therefore, The number of tapers used to form Pxx is 2*NW-1.
%   
%   Pxx = PMTM(X,NW,NFFT) specifies the FFT length used to calculate the
%   PSD estimates.  For real X, Pxx has (NFFT/2+1) rows if NFFT is even,
%   and (NFFT+1)/2 rows if NFFT is odd.  For complex X, Pxx always has
%   length NFFT.  If NFFT is specified as empty, NFFT is set to either
%   256 or the next power of 2 greater than the length of X, whichever is
%   larger.
%
%   [Pxx,W] = PMTM(X,NW,NFFT) returns the vector of normalized angular 
%   frequencies, W, at which the PSD is estimated.  W has units of 
%   radians/sample.  For real signals, W spans the interval [0,Pi] when
%   NFFT is even and [0,Pi) when NFFT is odd.  For complex signals, W 
%   always spans the interval [0,2*Pi).
%
%   [Pxx,W] = PMTM(X,NW,W) computes the two-sided PSD at the normalized
%   angular frequencies contained in vector W. W must have at least two
%   elements.
%
%   [Pxx,F] = PMTM(X,NW,NFFT,Fs) returns a PSD computed as a function of
%   physical frequency.  Fs is the sampling frequency specified in hertz.
%   If Fs is empty, it defaults to 1 Hz.
%
%   F is the vector of frequencies (in hertz) at which the PSD is
%   estimated.  For real signals, F spans the interval [0,Fs/2] when NFFT
%   is even and [0,Fs/2) when NFFT is odd.  For complex signals, F always
%   spans the interval [0,Fs).
%
%   [Pxx,F] = PMTM(X,NW,F,Fs) computes the two-sided PSD at the frequencies
%   contained in vector F.  F is a vector of frequencies in Hz with 2 or
%   more elements.
%
%   [Pxx,F] = PMTM(...,Fs,method) uses the algorithm specified in method 
%   for combining the individual spectral estimates:
%      'adapt'  - Thomson's adaptive non-linear combination (default).
%      'unity'  - linear combination with unity weights.
%      'eigen'  - linear combination with eigenvalue weights.
%
%   [Pxx,F] = PMTM(X,E,V,NFFT,Fs,method) is the PSD estimate, and frequency
%   vector from the data tapers in E and their concentrations V.  Type HELP
%   DPSS for a description of the matrix E and the vector V. By default,
%   PMTM drops the last eigenvector because its corresponding eigenvalue is
%   significantly smaller than 1.
%
%   [Pxx,F] = PMTM(X,DPSS_PARAMS,NFFT,Fs,method) uses the cell 
%   array DPSS_PARAMS containing the input arguments to DPSS (listed in
%   order, but excluding the first argument) to compute the data tapers. 
%   For example, PMTM(x,{3.5,'trace'},512,1000) calculates the prolate 
%   spheroidal sequences for NW=3.5, NFFT=512, and Fs=1000, and displays
%   the method that DPSS uses for this calculation. Type HELP DPSS for 
%   other options.
%
%   [Pxx,F] = PMTM(...,'DropLastTaper',DROPFLAG) specifies whether PMTM
%   should drop the last taper/eigenvector during the calculation. DROPFLAG
%   can be one of the following values: [ {true} | false ].
%       true  - the last taper/eigenvector is dropped 
%       false - the last taper/eigenvector is preserved
%
%   [Pxx,F,Pxxc] = PMTM(...,'ConfidenceLevel',P) returns the P*100%
%   confidence interval for Pxx, where P is a scalar between 0 and 1. The
%   default value for P is .95.  Confidence intervals are computed using a
%   chi-squared approach. Pxxc has twice as many columns as Pxx.
%   Odd-numbered columns contain the lower bounds of the confidence
%   intervals; even-numbered columns contain the upper bounds.  Thus,
%   Pxxc(M,2*N-1) is the lower bound and Pxxc(M,2*N) is the upper bound
%   corresponding to the estimate Pxx(M,N).
%
%   [...] = PMTM(X,...,FREQRANGE) returns the PSD over the specified range
%   of frequencies based upon the value of FREQRANGE:
%
%      'onesided' - returns the one-sided PSD of a real input signal X.
%         If NFFT is even, Pxx has length NFFT/2+1 and is computed over the
%         interval [0,pi].  If NFFT is odd, Pxx has length (NFFT+1)/2 and
%         is computed over the interval [0,pi). When Fs is specified, the
%         intervals become [0,Fs/2) and [0,Fs/2] for even and odd NFFT,
%         respectively.
%
%      'twosided' - returns the two-sided PSD for either real or complex
%         input X.  Pxx has length NFFT and is computed over the interval
%         [0,2*pi). When Fs is specified, the interval becomes [0,Fs).
%
%      'centered' - returns the centered two-sided PSD for either real or
%         complex X.  Pxx has length NFFT and is computed over the interval
%         (-pi, pi] for even length NFFT and (-pi, pi) for odd length NFFT.
%         When Fs is specified, the intervals become (-Fs/2, Fs/2] and
%         (-Fs/2, Fs/2) for even and odd length NFFT, respectively.
%
%      FREQRANGE may be placed in any position in the input argument list
%      after the second input argument, unless E and V are specified, in
%      which case FREQRANGE may be placed in any position after the third
%      input argument.  The default value of FREQRANGE is 'onesided' when X
%      is real and 'twosided' when X is complex.
%
%   PMTM(...) with no output arguments plots the PSD estimate (in decibels
%   per unit frequency) in the current figure window.
%
%   EXAMPLE:
%      Fs = 1000;   t = 0:1/Fs:.3;  
%      x = cos(2*pi*t*200)+randn(size(t)); % A cosine of 200Hz plus noise
%      pmtm(x,3.5,[],Fs);                  % Uses the default NFFT.
%
%   See also DPSS, PWELCH, PERIODOGRAM, PMUSIC, PBURG, PYULEAR, PCOV,
%   PMCOV, PEIG.

%   References: 
%     [1] Thomson, D.J."Spectrum estimation and harmonic analysis."
%         In Proceedings of the IEEE. Vol. 10 (1982). pp. 1055-1096.
%     [2] Percival, D.B. and Walden, A.T., "Spectral Analysis For Physical
%         Applications", Cambridge University Press, 1993, pp. 368-370. 

%   Copyright 1988-2014 The MathWorks, Inc.
%      

narginchk(1,13);
len = length(varargin);

% convert vectors to column vectors for backwards compatibility
if isvector(x)
  x = x(:);
end

if numel(size(x))>2
  error(message('signal:pmtm:MustBeVectorOr2DMatrix'));
end

tf = strcmpi('droplasttaper',varargin);
indx = find(tf==1);
if (~isempty(indx) && ~islogical(varargin{indx+1}))
    error(message('signal:pmtm:MustBeLogical', 'Droplasttaper'));
end

% If the 'droplasttaper' pv-pair is used, move it to the end of varargin
if (~isempty(indx) && (indx+1 ~= len))
    dummy = varargin(1:indx-1);
    dummy(indx:len-2) = varargin(indx+2:len);
    dummy(len-1:len) = varargin(indx:indx+1);
    varargin = dummy;
end

% Parse the inputs, set up default values, and return any error messages.
params = parseinputs(x,varargin{:});

% Cast to enforce Precision Rules
% Already checked for invalid character inputs (NW, NFFT,Fs) in 'parseinputs->psdoptions'
params.nfft = double(params.nfft);
params.Fs = double(params.Fs);

% Compute the two-sided power spectrum via MTM.
[S,k,w] = mtm_spectrum(x,params);

% Generate the freq vector in correct units to avoid roundoff errors due to
% converting frequency units later.
nfft = params.nfft;
[~,ncol] = size(nfft);

% Compute the 1-sided or 2-sided PSD [Power/freq] or mean-square [Power].
% Also, compute the corresponding freq vector & freq units.
[Pxx,f,units] = computepsd(S,w,params.range,params.nfft,params.Fs,'psd');  

% Calculate confidence limits ONLY when needed, since it can take a while.
% Enforce double precision arithmetic on the calculations. chi2conf already
% enforces double precision arithmetic.
if ~strcmp(params.conflevel,'omitted') && (nargout==0 || nargout>2)
  % 'ConfidenceLevel' pv-pair specified
  Nchan = size(x,2);
  c = chi2conf(params.conflevel,k);
  Pxxc(:,2:2:2*Nchan) = double(Pxx)*c(2);
  Pxxc(:,1:2:2*Nchan) = double(Pxx)*c(1);
elseif (~strcmp(params.ConfInt,'omitted') && nargout==0) || nargout>2
  % using legacy syntax
  if ~strcmp(params.ConfInt,'omitted')
    c = chi2conf(params.ConfInt,k);
  else % use default value of .95
    c = chi2conf(.95,k);
  end      
  Nchan = size(x,2);
  Pxxc(:,2:2:2*Nchan) = double(Pxx)*c(2);
  Pxxc(:,1:2:2*Nchan) = double(Pxx)*c(1);
else
  Pxxc = [];
end

% Output
if nargout==0
   % If no output arguments are specified plot the PSD w/ conf intervals.
   f = {f};
   if strcmpi(units,'Hz'), f = [f {'Fs',params.Fs}]; end
   hpsd = dspdata.psd([Pxx Pxxc],f{:},'SpectrumType',params.range);

   % Create a spectrum object to store in the PSD object's metadata.
   hspec = spectrum.mtm(params.E,params.V,params.MTMethod);
   hpsd.Metadata.setsourcespectrum(hspec);
      
   if params.centerdc
     centerdc(hpsd);
   end

   hp = plot(hpsd);
   if 3*size(x,2)==numel(hp)
     nChan = size(x,2);
     color = get(hp,'Color');
     for i=1:nChan
       set(hp(nChan+2*i-1), 'Color',color{i}, 'LineStyle','-.');
       set(hp(nChan+2*i), 'Color',color{i}, 'LineStyle','-.');
     end
   end     
   return
end

% center dc if needed
if params.centerdc
   [Pxx, f, Pxxc] = psdcenterdc(Pxx, f, Pxxc, params);
end

% If the input is a vector and a row frequency vector was specified,
% return output as a row vector for backwards compatibility.
if ncol > 1 && nargout > 0 && isvector(x)
    Pxx = Pxx.';
    f = f.';
    % preserve (undocumented) behavior with legacy syntax.
    if strcmp(params.conflevel,'omitted') && nargout >= 3
        Pxxc = Pxxc.';
    end
end

if isa(Pxx,'single')
  % Cast to enforce precision rules.
  f = single(f);
end  

if nargout==1
   varargout = {Pxx};
elseif nargout==2
   varargout = {Pxx,f};
elseif nargout==3
   if ~strcmp(params.conflevel,'omitted')
      % use preferred output order
      varargout = {Pxx,f,Pxxc};
   else
      % use legacy output order
      varargout = {Pxx,Pxxc,f};
   end
end

%----------------------------------------------------------------------
function [S,k,w] = mtm_spectrum(x,params)
%MTM_SPECTRUM Compute the power spectrum via MTM.
%
% Inputs:
%   x      - Input data vector.
%   params - Structure containing pmtm's input parameter list, except for
%            the input data sequence, x; it contains the following fields:
%      nfft     - Number of frequency points to evaluate the PSD at; 
%                 the default is max(256,2^nextpow2(N)).
%      Fs       - The sampling frequency; default is 1.
%      range    - default is 'onesided' or real signals and 'twosided' for 
%               - complex signals.
%      ConfInt  - Confidence interval; default is .95.
%      MTMethod - Algorithm used in MTM; default is 'adapt'.
%      E        - Matrix containing the discrete prolate spheroidal 
%                 sequences (dpss).
%      V        - Vector containing the concentration of the dpss.
%      NW       - Time-bandwidth product; default is 4.
%
% Outputs:
%   S      - Power spectrum computed via MTM.
%   k      - Number of sequences used to form Pxx.
%   w      - Frequency vector for which DFT is calculated

% Extract some parameters from the input structure for convenience.
nfft = params.nfft;
E  = params.E;
V  = params.V;
Fs = params.Fs;

if isempty(Fs)
  Fs = 2*pi;
end

N = size(x,1);
Nchan = size(x,2);
k = length(V);

if length(nfft) > 1, 
    isfreqVector = true;     
    nfft_mod = length(nfft);
else 
    isfreqVector = false;
    nfft_mod = nfft;
end

if isa(x,'single') || isa(E,'single')
  S = zeros(nfft_mod, Nchan, 'single');
else
  S = zeros(nfft_mod, Nchan);
end

for chan=1:Nchan
    % Compute the windowed DFTs.
    if (~isfreqVector && N<=nfft) || isfreqVector 

        % Compute DFT using FFT or Goertzel
        [Xx,w] = computeDFT(bsxfun(@times,E(:,1:k),x(:,chan)),nfft(:),Fs);    
        Sk = abs(Xx).^2;

    else % Wrap the data modulo nfft if N > nfft. Note we cannot use datawrap 
        % and FFT because datawrap does not support matrices
        % use CZT to compute DFT on nfft evenly spaced samples around the
        % unit circle:
        Sk = abs(czt(bsxfun(@times,E(:,1:k),x(:,chan)),nfft(:))).^2;
        w = psdfreqvec('npts',nfft,'Fs',Fs);
    end

    % Compute the MTM spectral estimates, compute the whole spectrum 0:nfft.
    switch params.MTMethod,

    case 'adapt'
       % Set up the iteration to determine the adaptive weights: 

       sig2=x(:,chan)'*x(:,chan)/N;              % Power
       Schan=(Sk(:,1)+Sk(:,2))/2;    % Initial spectrum estimate   
       S1=zeros(nfft_mod,1);  

       % The algorithm converges so fast that results are
       % usually 'indistinguishable' after about three iterations.

       % This version uses the equations from [2] (P&W pp 368-370).

       % Set tolerance for acceptance of spectral estimate:
       tol=.0005*sig2/nfft_mod;
       a=bsxfun(@times,sig2,(1-V));

       % Do the iteration:
       while sum(abs(Schan-S1)/nfft_mod)>tol
          % calculate weights
          b=(Schan*ones(1,k))./(Schan*V'+ones(nfft_mod,1)*a'); 
          % calculate new spectral estimate
          wk=(b.^2).*(ones(nfft_mod,1)*V');
          S1=sum(wk'.*Sk')./ sum(wk,2)';
          S1=S1';
          Stemp=S1; S1=Schan; Schan=Stemp;  % swap S and S1
       end
    case {'unity','eigen'}
       % Compute the averaged estimate: simple arithmetic averaging is used. 
       % The Sk can also be weighted by the eigenvalues, as in Park et al. 
       % Eqn. 9.; note that the eqn. apparently has a typo; as the weights
       % should be V and not 1/V.
       if strcmp(params.MTMethod,'eigen')
          wt = V(:);    % Park estimate
       else
          wt = ones(k,1);
       end
       Schan = Sk*wt/k;
    end
    S(:,chan) = Schan;
end

%----------------------------------------------------------------------
function params = parseinputs(x,varargin)
%PARSEINPUTS Parse the inputs passed to pmtm.m and return a structure
%            containing all the parameters passed to PMTM set to either
%            default values or user defined values.
%
% Inputs:
%   x        - Input data vector.
%   varargin - Input parameter list passed to pmtm, except for x.
%
% Outputs:
%   params   - Structure containing pmtm's input parameter list, except for
%              the input data sequence, x; it contains the following fields:
%      nfft     - Number of frequency points to evaluate the PSD at; 
%                 the default is max(256,2^nextpow2(N)).
%      Fs       - The sampling frequency; default is .
%      range    - default is 'onesided' or real signals and 'twosided' for 
%               - complex signals.
%      conflevel- Confidence level (preferred syntax)
%      ConfInt  - Confidence interval; default is .95. (legacy syntax)
%      MTMethod - Algorithm used in MTM; default is 'adapt'.
%      E        - Matrix containing the discrete prolate spheroidal 
%                 sequences.
%      V        - Vector containing the concentration of the dpss.
%      NW       - Time-bandwidth product; default is 4.
%
%   err_msg  - String containing an error message if an error occurred.

if any(strcmp(varargin, 'whole'))
    warning(message('signal:pmtm:invalidRange', 'whole', 'twosided'));
elseif any(strcmp(varargin, 'half'))
    warning(message('signal:pmtm:invalidRange', 'half', 'onesided'));
end

% Set default parameter values.
N = size(x,1);
params  = [];

% Parse the input arguments up to NFFT (exclusive). 
% If E and V are not specified, calculate them.
[E,V,NW,indx,nfft_temp,varargin] = getEV(N,varargin{:});

% Cast to enforce Precision Rules
if (any([signal.internal.sigcheckfloattype(x,'single','pmtm','X') ...
    signal.internal.sigcheckfloattype(E,'single','pmtm','E') ...
    signal.internal.sigcheckfloattype(V,'single','pmtm','V')]))
  x = single(x);
  E = single(E);
  V = single(V);
end
NW = double(NW);

if isreal(x) && (length(nfft_temp) <= 1), 
   range = 'onesided';
else
   range = 'twosided'; 
end

% NOTE: The psdoptions function REQUIRES a structure with the following 
%       fields.  Any changes to the structure, such as adding/removing 
%       fields, should be done after the call to psdoptions.
params.nfft    = max(256,2^nextpow2(N));
params.Fs      = [];
params.range   = range;
params.centerdc = false;
params.conflevel = 'omitted';

params.ConfInt = 'omitted';
params.MTMethod= 'adapt';

% Call psdoptions to handle the remaining input arg list starting with NFFT.
% Overwrite default options with user specified options (if specified).
if indx <= numel(varargin)
  % Invalid character inputs for NW, NFFT, W, E,V and Fs is checked here
   [params,err_msg,err_msgobj] = psdoptions(isreal(x),params,varargin{indx:end});
   if err_msg, error(err_msgobj), end;     
   
   if ~strcmp(params.conflevel,'omitted') && ~strcmp(params.ConfInt, 'omitted')
       % user specified too many scalar inputs in conjunction with 'ConfidenceLevel'
       error(message('signal:pmtm:TooManyScalarNumericInputs'));
   end
   
   if length(params.nfft) > 1,
       if strcmpi(params.range,'onesided')
           warning(message('signal:pmtm:InconsistentRangeOption'));
       end
       params.range = 'twosided';
   end
end

% Add remaining fields to the return structure.
params.E  = E;
params.V  = V;
params.NW = NW;

%----------------------------------------------------------------------
function [E,V,NW,indx,nfft_temp,varargin] = getEV(N,varargin)
% GETEV  Parse the input arguments up to, but not including, Nfft and 
%        calculate E and V if not specified.
%
% Inputs:
%   N        - Length of the input data sequence, x.
%   varargin - Input parameter list passed to pmtm, except for x.
%
% Outputs:
%   E        - Matrix containing the discrete prolate spheroidal 
%              sequences (dpss).
%   V        - Vector containing the concentration of the dpss.
%   NW       - Time-bandwidth product; default is 4.
%   indx     - Index indicating starting location of options in pmtm's 
%              input argument list.
%   nfft_temp - NFFT or Frequency vector specified. Empty if not specified 

% Define defaults & initialize output variables (in case of early return).
NW      = 4;
indx    = 1;  % Index where the options begin in the input arg list
nfft_temp = [];

tf = strcmpi('droplasttaper',varargin);
loc = find(tf==1);
if ~isempty(loc)
    dlt = varargin{loc+1};     % droplasttaper
    varargin = varargin(1:loc-1);
else
    dlt = true;              % default value
end

% The 2nd input arg to pmtm can be a
%    1. (X,NW,...)            scalar
%    2. (X,E,V,...)           matrix E, hence, 3rd input must be a vector (V) 
%    3. (X,{dpss_params},...) cell containing the input argument list to dpss 
%    4. a string specifying a psdoption
if ~isempty(varargin) && ~ischar(varargin{1})
   if ~isempty(varargin{1})
       NW = varargin{1};
   end
   indx = 2;
   if iscell(NW),           % NW is a cell array => dpss_params
      if (NW{1}<1.25 && dlt)
          error(message('signal:pmtm:insufficientTimebandwidthproduct', 'NW', '1.25', 'Droplasttaper', 'true'));
      end 
      if (NW{1}<0.75 && ~dlt)
          error(message('signal:pmtm:insufficientTimebandwidthproduct', 'NW', '0.75', 'Droplasttaper', 'false'));          
      end       
      [E,V] = dpss(N,NW{:}); 
      numvec = length(V);
      if dlt
           if numvec > 2
               E = E(:,1:numvec-1);
               V = V(1:numvec-1);
           else
               error(message('signal:pmtm:inadequateNumtapers', '3', 'Droplasttaper', 'true'));
           end
      else
           if numvec < 2
               error(message('signal:pmtm:inadequateNumtapers', '2', 'Droplasttaper', 'false'));
           end
      end
      NW = NW{1};
      if nargin > 2, nfft_temp = findNFFT(varargin{2:end}); end
   elseif length(NW)>1,     % NW is the matrix E (==>V must be specified)
      E = NW;
      if length(varargin)<2,
         error(message('signal:pmtm:MustProvideVWithE', 'V', 'E'));
      else
         V = varargin{2};
         if nargin > 3, nfft_temp = findNFFT(varargin{3:end}); end
      end
      numvec = length(V);
      if size(E,2)~=numvec
         error(message('signal:pmtm:MismatchEV', 'E', 'V'));
      end     
      NW = size(E,2)/2;  
      indx = 3; % Update index of beginning of options in the input arg list      
      if dlt
          if (numvec < 3)
               error(message('signal:pmtm:inadequateNumtapers', '3', 'Droplasttaper', 'true'));
          else
              E = E(:,1:numvec-1);
              V = V(1:numvec-1);
          end
      else
          if(numvec < 2)
               error(message('signal:pmtm:inadequateNumtapers', '2', 'Droplasttaper', 'false'));
          end
      end
   else                      % NW is a scalar
       if (NW<1.25 && dlt)
          error(message('signal:pmtm:insufficientTimebandwidthproduct', 'NW', '1.25', 'Droplasttaper', 'true'));
       end
       if (NW<0.75 && ~dlt)
          error(message('signal:pmtm:insufficientTimebandwidthproduct', 'NW', '0.75', 'Droplasttaper', 'false'));
       end
       % Get the dpss, one way or another:
       [E,V] = dpss(N,NW);
       numvec = length(V);
       if dlt
           if numvec > 2
               E = E(:,1:numvec-1);
               V = V(1:numvec-1);
           else
               error(message('signal:pmtm:inadequateNumtapers', '3', 'Droplasttaper', 'true'));
           end
       else
           if numvec < 2
               error(message('signal:pmtm:inadequateNumtapers', '2', 'Droplasttaper', 'false'));
           end
       end
       if nargin > 2, nfft_temp = findNFFT(varargin{2:end}); end
   end
else
   % Get the dpss, one way or another:
   [E,V] = dpss(N,NW);
   numvec = length(V);
      if dlt
           if numvec > 2
               E = E(:,1:numvec-1);
               V = V(1:numvec-1);
           else
               error(message('signal:pmtm:inadequateNumtapers', '3', 'Droplasttaper', 'true'));
           end
      else
           if numvec < 2
               error(message('signal:pmtm:inadequateNumtapers', '2', 'Droplasttaper', 'false'));
           end
      end
   nfft_temp = [];
end


%------------------------------------------------------------------
function nfftTemp = findNFFT(varargin)
% FINDNFFT Finds the specified NFFT or frequency vector from the optional
% arguments passed

nfftTemp = [];
for cnt = 1:length(varargin)
    if isnumeric(varargin{cnt}), 
        nfftTemp = varargin{cnt};
        break;
    end
end

% [EOF] pmtm.m




