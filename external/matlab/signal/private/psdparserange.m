function [Pxx, F, Frange, rbw, extraArgs, status] = psdparserange(funcName, beta, varargin)
%PSDMASKACTIVERANGE Helper function for ranged power and psd estimates
%   computes and/or returns a PSD for use based upon the input
%
%      funcName - name of calling function to use with error messages
%
%      beta - coefficient for Kaiser window.
%             Usually 0 (rectangular) or 38 (for use with sinusoids)
%   
%      varargin - supports:
%                    FUNC(X) 
%                    FUNC(X, Fs)
%                    FUNC(Pxx, F)
%                    FUNC(Sxx, F, RBW)
%                    FUNC(..., FREQRANGE)
%                    FUNC(..., FREQRANGE, EXTRAARGS)
%          X  - time vector intput
%          Fs - sample rate
%          F - frequency vector
%          Pxx - PSD estimate
%          Sxx - spectral power estimate
%          FREQRANGE - must be empty or a two-element vector
%          EXTRAARGS - anything following the FREQRANGE
%          NORMF - indicates if plot is logical.
%   Copyright 2014 The MathWorks, Inc.

n = numel(varargin);

oneSided = false;
hasNyquist = false;

if n<2 || isempty(varargin{2}) || isscalar(varargin{2})
  inputType = 'time';
  oneSided = true;
  [Pxx, F, Frange, rbw] = parseTime(funcName, beta, varargin{1:min(n,3)});
  hasNyquist = hasNyquistBin(varargin{1});
  extraArgs = varargin(4:end);
elseif n<3 || ~isscalar(varargin{3})
  inputType = 'psd';
  [Pxx, F, Frange, rbw] = parsePSD(funcName, varargin{1:min(n,3)});
  extraArgs = varargin(4:end);
else
  inputType = 'power';
  [Pxx, F, Frange, rbw] = parsePower(funcName, varargin{1:min(n,4)});
  extraArgs = varargin(5:end);
end

% report normalized frequency flag
status.normF = numel(varargin)==1 || isequal(varargin{2},[]);
status.inputType = inputType;
status.oneSided = oneSided;
status.hasNyquist = hasNyquist;

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function flag = hasNyquistBin(x)
if isvector(x)
  x = x(:);
end
n = size(x,1);
flag = n/2 == round(n/2);

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function [Pxx, F, Frange, rbw] = parseTime(funcName, beta, x, fs, Frange)
% force column vector before checking attributes
if isvector(x)
  x = x(:);
end

% validate x
validateattributes(x,{'numeric'},{'2d','finite'}, funcName,'x',1);

% validate fs
if nargin > 3 && ~isempty(fs)
  validateattributes(fs, {'numeric'},{'real','finite','scalar','positive'}, ...
    funcName,'Fs',2);
else
  % use angular frequency
  fs = 2*pi;
end

% use Kaiser window to reduce effects of leakage
n = size(x,1);
w = kaiser(n,beta);
rbw = enbw(w,fs);

% use one-sided PSD for real signals, otherwise use centered for complex
if isreal(x)
  [Pxx, F] = periodogram(x,w,n,fs,'psd');
else
  [Pxx, F] = periodogram(x,w,n,fs,'centered','psd');
end

% check freq range vector.  If not specified, return empty.
if nargin > 4 && ~isempty(Frange)
  if fs==2*pi
    validateattributes(Frange, {'numeric'},{'real','finite','increasing','size',[1 2]}, ...
      funcName,'FREQRANGE',2);
  else
    validateattributes(Frange, {'numeric'},{'real','finite','increasing','size',[1 2]}, ...
      funcName,'FREQRANGE',3);
  end
  if Frange(1)<F(1) || Frange(2)>F(end)
    error(message('signal:psdparserange:FreqRangeOutOfBounds'));
  end
else
  Frange = [];
end

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function [Pxx, F, Frange, rbw] = parsePSD(funcName, Pxx, F, Frange)
% force column vector before checking attributes
if isvector(Pxx)
  Pxx = Pxx(:);
end

if isvector(F)
  F = F(:);
end

% validate Pxx
validateattributes(Pxx,{'numeric'},{'2d','nonnegative'}, funcName,'Pxx',1);

% validate F
validateattributes(F,{'numeric'},{'real','vector','finite'}, funcName,'F',2);

if size(Pxx,1) ~= numel(F)
    error(message('signal:psdparserange:FreqVectorMismatch'));
end

% check freq range vector.  If not specified, return empty.
if nargin > 3 && ~isempty(Frange)
  validateattributes(Frange, {'numeric'},{'real','finite','increasing','size',[1 2]}, ...
    funcName,'FREQRANGE',3);
  if Frange(1)<F(1) || Frange(2)>F(end)
    error(message('signal:psdparserange:FreqRangeOutOfBounds'));
  end
else
  Frange = [];
end

% no need to return an RBW for PSD.
rbw = NaN;

function [Pxx,F,Frange,rbw] = parsePower(funcName, Sxx, F, rbw, Frange)

% force column vector before checking attributes
if isvector(Sxx)
  Sxx = Sxx(:);
end

if isvector(F)
  F = F(:);
end

% validate Sxx
validateattributes(Sxx,{'numeric'},{'2d','nonnegative'}, funcName,'Sxx',1);

% validate F
validateattributes(F,{'numeric'},{'real','vector','finite'}, funcName,'F',2);
if size(Sxx,1) ~= numel(F)
    error(message('signal:psdparserange:FreqVectorMismatch'));
end

% ensure specified RBW is larger than a bin width
df = mean(diff(F));

validateattributes(rbw,{'double'},{'real','finite','positive','scalar','>=',df}, ...
    funcName,'RBW',3);
  
% check freq range vector.  If not specified, return empty.
if nargin > 4 && ~isempty(Frange)
  validateattributes(Frange, {'numeric'},{'real','finite','increasing','size',[1 2]}, ...
    funcName,'FREQRANGE',3);
  if Frange(1)<F(1) || Frange(2)>F(end)
    error(message('signal:psdparserange:FreqRangeOutOfBounds'));
  end
else
  Frange = [];
end
Pxx = Sxx/rbw;

  
