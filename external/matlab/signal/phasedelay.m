function varargout = phasedelay(b,varargin)
%PHASEDELAY Phase delay of digital filter
%   [PHI,W] = PHASEDELAY(B,A,N) returns the N-point phase delay response
%   vector PHI (in samples) and the N-point frequency vector W (in
%   radians/sample) of the filter:
%
%               jw               -jw              -jmw
%        jw  B(e)    b(1) + b(2)e + .... + b(m+1)e
%     H(e) = ---- = ------------------------------------
%               jw               -jw              -jnw
%            A(e)    a(1) + a(2)e + .... + a(n+1)e
%
%   given numerator and denominator coefficients in vectors B and A.
%
%   [PHI,W] = PHASEDELAY(SOS,N) computes the phase delay response of the
%   filter specified using the second order sections matrix SOS. SOS is a
%   Kx6 matrix, where the number of sections, K, must be greater than or
%   equal to 2. Each row of SOS corresponds to the coefficients of a second
%   order filter. From the transfer function displayed above, the ith row
%   of the SOS matrix corresponds to [bi(1) bi(2) bi(3) ai(1) ai(2) ai(3)].
%
%   [PHI,W] = PHASEDELAY(D,N) computes the phase delay response of the
%   digital filter, D. You design a digital filter, D, by calling the 
%   <a href="matlab:help designfilt">designfilt</a> function.
%
%   In all cases, the phase response is evaluated at N points equally
%   spaced around the upper half of the unit circle. If N isn't specified,
%   it defaults to 512.
%
%   [PHI,W] = PHASEDELAY(...,N,'whole') uses N points around the whole unit
%   circle.
%
%   [PHI,F] = PHASEDELAY(...,N,Fs) and [PHI,F] =PHASEDELAY(...,N,'whole',Fs)
%   return a frequency vector, F, in Hz when you specify the sample rate Fs
%   in Hz.
%
%   PHI = PHASEDELAY(...,W) and PHI = PHASEDELAY(..,F,Fs) return the phase
%   delay response evaluated at the points specified in frequency vectors W
%   (in radians/sample), or F (in Hz).
%
%   PHASEDELAY(...) with no output arguments plots the phase delay response
%   of the filter in the current figure window.
%
%   % Example 1:
%   %   Design a lowpass FIR filter with normalized cut-off frequency at
%   %   0.3 and determine its phase delay.
%
%   b=fircls1(54,0.3,0.02,0.008);
%   phasedelay(b)
%
%   % Example 2:
%   %   Design a 5th order lowpass elliptic IIR filter and determine its
%   %   phase delay.
%
%   [b,a] = ellip(5,0.5,20,0.4);
%   phasedelay(b,a)
%
%   % Example 3:
%   %   Design a Butterworth highpass IIR filter, represent its coefficients
%   %   using second order sections, and display its phase delay response.
%
%   [z,p,k] = butter(6,0.7,'high');
%   SOS = zp2sos(z,p,k);
%   phasedelay(SOS)
%
%   % Example 4:
%   %   Use the designfilt function to design a highpass IIR digital filter 
%   %   with order 8, passband frequency of 75 KHz, and a passband ripple 
%   %   of 0.2 dB. Sample rate is 200 KHz. Visualize the phase delay response 
%   %   using 2048 frequency points.
%  
%   D = designfilt('highpassiir', 'FilterOrder', 8, ...
%            'PassbandFrequency', 75e3, 'PassbandRipple', 0.2,...
%            'SampleRate', 200e3);
%
%   phasedelay(D,2048)
%
%   See also FREQZ, PHASEZ, ZEROPHASE, GRPDELAY and FVTOOL.

%   Copyright 1988-2013 The MathWorks, Inc.

narginchk(1,5);

isTF = true; % True if dealing with a transfer function

if all(size(b)>[1 1])
    % Input is a matrix, check if it is a valid SOS matrix
    if size(b,2) ~= 6
        error(message('signal:signalanalysisbase:invalidinputsosmatrix'));
    end
    % Checks if SOS is a valid numeric data input
    signal.internal.sigcheckfloattype(b,'','phasedelay','SOS');
    
    isTF = false; % SOS instead of transfer function
end

if isTF
    if isempty(varargin)
        a = 1; % Assume FIR
    else
        a = varargin{1};
        varargin(1) = [];
    end
    
    % Checks if B and A are valid numeric data inputs
    signal.internal.sigcheckfloattype(b,'','phasedelay','B');
    signal.internal.sigcheckfloattype(a,'','phasedelay','A');

    if ~isreal(b) || ~isreal(a),
        [phi,w,s] = phasez(b,a,varargin{:});
    else
        % Use the continuous phase
        [~,w,phi,s] = zerophase(b,a,varargin{:});
    end
else
    [~,w, phi, s] = zerophase(b, varargin{:});
end

% Note that phi and w span between [0, pi)/[0, fs/2)
phd = dividenowarn(-phi,w);

% Cast to enforce precision rules
% Only cast if output is requested, otherwise, plot using double precision
% frequency vector. 
if nargout > 1 && isa(phd,'single')
  w = single(w);
end

% Parse outputs
switch nargout
    case 0
        % Plot when no output arguments are given
        phasedelayplot(phd,w,s);
    case 1
        varargout = {phd};
    case 2
        varargout = {phd,w};
    case 3
        varargout = {phd,w,s};
end


%-------------------------------------------------------------------------------
function phasedelayplot(phd,w,s)

% Cell array of the standard frequency units strings (used for the Xlabels)
frequnitstrs = getfrequnitstrs;
if isempty(s.Fs),
    xlab = frequnitstrs{1};
    w    = w./pi; % Scale by pi in the plot
    ylab = 'Phase delay (samples)';
else
    xlab = frequnitstrs{2};
    ylab = 'Phase delay (rad/Hz)';
end

plot(w,phd);
xlabel(xlab);
ylabel(ylab);
grid on;

% [EOF]
