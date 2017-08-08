function varargout = phasez(b,varargin)
%PHASEZ Phase response of digital filter
%   [PHI,W] = PHASEZ(B,A,N) returns the N-point unwrapped phase response
%   vector PHI and the N-point frequency vector W in radians/sample of
%   the filter:
%               jw               -jw              -jmw
%        jw  B(e)    b(1) + b(2)e + .... + b(m+1)e
%     H(e) = ---- = ------------------------------------
%               jw               -jw              -jnw
%            A(e)    a(1) + a(2)e + .... + a(n+1)e
%   given numerator and denominator coefficients in vectors B and A. 
%
%   [PHI,W] = PHASEZ(SOS,N) returns the N-point unwrapped phase response
%   given the second order sections matrix SOS. SOS is a Kx6 matrix, where
%   the number of sections, K, must be greater than or equal to 2. Each row
%   of SOS corresponds to the coefficients of a second order filter. From
%   the transfer function displayed above, the ith row of the SOS matrix
%   corresponds to [bi(1) bi(2) bi(3) ai(1) ai(2) ai(3)].
%
%   [PHI,W] = PHASEZ(D,N) returns the N-point unwrapped phase response
%   given the digital filter, D. You design a digital filter, D, by calling
%   the <a href="matlab:help designfilt">designfilt</a> function.
%
%   In all cases, the phase response is evaluated at N points equally
%   spaced around the upper half of the unit circle. If N isn't specified,
%   it defaults to 512.
%
%   [PHI,W] = PHASEZ(...,N,'whole') uses N points around the whole unit
%   circle.
%
%   PHI = PHASEZ(...,W) returns the phase response at frequencies
%   designated in vector W, in radians/sample (normally between 0 and pi).
%
%   [PHI,F] = PHASEZ(...,N,Fs) and [PHI,F] = PHASEZ(...,N,'whole',Fs)
%   return phase vector F (in Hz), where Fs is the sampling frequency (in
%   Hz).
%
%   PHI = PHASEZ(...,F,Fs) returns the phase response at the frequencies
%   designated in vector F (in Hz), where Fs is the sampling frequency (in
%   Hz).
%
%   PHASEZ(...) with no output arguments plots the unwrapped phase of
%   the filter.
%
%   % Example 1:
%   %   Design a lowpass FIR filter with normalized cut-off frequency at 
%   %   0.3 and determine its phase response.
%
%   b=fircls1(54,0.3,0.02,0.008);
%   phasez(b)                       
%
%   % Example 2: 
%   %   Design a 5th order lowpass elliptic IIR filter and determine its
%   %   phase response.
%
%   [b,a] = ellip(5,0.5,20,0.4);
%   phasez(b,a,512,'whole');        
%
%   % Example 3:
%   %   Design a Butterworth highpass IIR filter, represent its coefficients
%   %   using second order sections, and display its phase response.
%
%   [z,p,k] = butter(6,0.7,'high');
%   SOS = zp2sos(z,p,k);    
%   phasez(SOS)      
%
%   % Example 4:
%   %   Use the designfilt function to design a highpass IIR digital filter 
%   %   with order 8, passband frequency of 75 KHz, and a passband ripple 
%   %   of 0.2 dB. Sample rate is 200 KHz. Visualize the phase response 
%   %   using 2048 frequency points.
%  
%   D = designfilt('highpassiir', 'FilterOrder', 8, ...
%            'PassbandFrequency', 75e3, 'PassbandRipple', 0.2,...
%            'SampleRate', 200e3);
%
%   phasez(D,2048)
%
%   See also FREQZ, PHASEDELAY, GRPDELAY and FVTOOL.

%   Copyright 1988-2013 The MathWorks, Inc.

narginchk(1,5);

a = 1; % Assume FIR for now
isTF = true; % True if dealing with a transfer function

if all(size(b)>[1 1])
  % Input is a matrix, check if it is a valid SOS matrix
  if size(b,2) ~= 6
    error(message('signal:signalanalysisbase:invalidinputsosmatrix'));
  end
  % Checks if SOS is a valid numeric data input
  signal.internal.sigcheckfloattype(b,'','phasez','SOS');
  isTF = false; % SOS instead of transfer function  
end

if isTF
  if nargin > 1
    a = varargin{1};
    varargin(1) = [];
  end

  if all(size(a)>[1 1])
    error(message('signal:signalanalysisbase:inputnotsupported'));
  end
end

% Define the new N-point frequency vector where the frequency response is evaluated
[upn_or_w, upfactor, iswholerange, options, addpoint] = findfreqvector(varargin{:});

% Compute the frequency response (freqz)
if isTF
  % Checks if B and A are valid numeric data inputs
  signal.internal.sigcheckfloattype(b,'','phasez','B');
  signal.internal.sigcheckfloattype(a,'','phasez','A');
   
  % freqz casts inputs without single precision bearing (N, F and Fs) to
  % double
  [h, w, s, options] = freqz(b, a, upn_or_w, options{:});
  
  % Extract phase from frequency response
  [phi,w] = extract_phase(h,upn_or_w,iswholerange,upfactor,w,addpoint);
else
  %Compute phase for individual sections, not on the entire SOS matrix to
  %avoid close-to-zero responses. Add phase of individual sections.  
  [h, w, s, options1] = freqz(b(1,1:3), b(1,4:6), upn_or_w, options{:});
  [phi,w] = extract_phase(h,upn_or_w,iswholerange,upfactor,w,addpoint);
  for indx = 2:size(b, 1)
    h = freqz(b(indx,1:3), b(indx,4:6), upn_or_w, options{:});
    phi = phi + extract_phase(h,upn_or_w,iswholerange,upfactor,w,addpoint);            
  end
  options = options1;
end

% Update options
options.nfft = length(w);
options.w    = w;

if length(a)==1 && length(b)==1,
    % Scalar case
    phi = angle(b/a)*ones(length(phi),1);
end

% Cast to enforce precision rules
% Only cast if output is requested, otherwise, plot using double precision
% frequency vector. 
if nargout > 1 && isa(phi,'single')
  w = single(w);
end

% Parse outputs
switch nargout
    case 0
        % Plot when no output arguments are given
        phaseplot(phi,w,s);
    case 1
        varargout = {phi};
    case 2
        varargout = {phi,w};
    case 3
        varargout = {phi,w,s};
    case 4
        varargout = {phi,w,s,options};
end

%-------------------------------------------------------------------------------
function phaseplot(phi,w,s)

% Cell array of the standard frequency units strings (used for the Xlabels)
frequnitstrs = getfrequnitstrs;
switch lower(s.xunits),
    case 'rad/sample',
        xlab = frequnitstrs{1};
        w    = w./pi; % Scale by pi in the plot
    case 'hz',
        xlab = frequnitstrs{2};
    case 'khz',
        xlab = frequnitstrs{3};
    case 'mhz',
        xlab = frequnitstrs{4};
    case 'ghz',
        xlab = frequnitstrs{5};
    otherwise
        xlab = s.xunits;
end

plot(w,phi);
title(getString(message('signal:phasez:PhaseResponse')));
xlabel(xlab);
ylabel(getString(message('signal:phasez:Phaseradians')));
grid on;

% [EOF]
