function [P,f] = computeperiodogram(x,win,nfft,esttype,varargin)
%COMPUTEPERIODOGRAM   Periodogram spectral estimation.
%   This function is used to calculate the Power Spectrum Sxx, and the
%   Cross Power Spectrum Sxy.
%
%   Sxx = COMPUTEPERIODOGRAM(X,WIN,NFFT) where x is a vector returns the
%   Power Spectrum over the whole Nyquist interval, [0, 2pi).
%
%   Sxy = COMPUTEPERIODOGRAM({X,Y},WIN,NFFT) returns the Cross Power
%   Spectrum over the whole Nyquist interval, [0, 2pi).
%
%   Inputs:
%    X           - Signal vector or a cell array of two elements containing
%                  two signal vectors.
%    WIN         - Window
%    NFFT        - Number of frequency points (FFT) or vector of
%    frequencies at which periodogram is desired (Goertzel)
%    WINCOMPFLAG - A string indicating the type of window compensation to
%                  be done. The choices are: 
%                  'ms'    - compensate for Mean-square (Power) Spectrum;
%                            maintain the correct power peak heights.
%                  'power' - compensate for Mean-square (Power) Spectrum;
%                            maintain the correct power peak heights.
%                  'psd'   - compensate for Power Spectral Density (PSD);
%                            maintain correct area under the PSD curve.
%
%   Output:
%    Sxx         - Power spectrum [Power] over the whole Nyquist interval. 
%      or
%    Sxy         - Cross power spectrum [Power] over the whole Nyquist
%                  interval.

%   Author(s): P. Pacheco
%   Copyright 1988-2014 The MathWorks, Inc.

narginchk(3,5);
if nargin < 4,
    esttype = 'psd'; % Default, compensate for window's power.
end

if nargin < 5 || isempty(varargin{1}),
    Fs = 2*pi;
else
    Fs = varargin{1};     
end

% Validate inputs and convert row vectors to column vectors.
[x,~,y,is2sig,win] = validateinputs(x,win,nfft);

% Window the data
xw = bsxfun(@times,x,win);
if is2sig, yw = bsxfun(@times,y,win); end 

% Evaluate the window normalization constant.  A 1/N factor has been
% omitted since it will cancel below.
if any(strcmpi(esttype,{'ms','power'}))
    % The window is convolved with every power spectrum peak, therefore
    % compensate for the DC value squared to obtain correct peak heights.
    U = sum(win)^2;
else
    U = win'*win;  % compensates for the power of the window.
end

% Compute the periodogram power spectrum [Power] estimate
% A 1/N factor has been omitted since it cancels

[Xx,f] = computeDFT(xw,nfft,Fs);
if is2sig, [Yy,f] = computeDFT(yw,nfft,Fs); end

P = Xx.*conj(Xx)/U;      % Auto spectrum.
if is2sig,
    P = bsxfun(@times,Xx,conj(Yy))/U;  % Cross spectrum.
end

%--------------------------------------------------------------------------
function [x,Lx,y,is2sig,win] = validateinputs(x,win,~)
% Validate the inputs to computexperiodogram and convert row vectors to
% column vectors for backwards compatiblility with R2014a and prior
% releases

% Set defaults and convert to row vectors to columns.
y     = [];
is2sig= false;
win   = win(:);
Lw    = length(win);

% Determine if one or two signal vectors was specified.
if iscell(x),
    if length(x) > 1,
        y = x{2};
        if isvector(y)
            y = y(:);
        end
        is2sig = true;
    end
    x = x{1};
end

if isvector(x)
    x = x(:);
end

Lx = size(x,1);

if is2sig,
    Ly  = size(y,1);
    if Lx ~= Ly,
        error(message('signal:computeperiodogram:invalidInputSignalLength'))
    end
    if size(x,2)~=1 && size(y,2)~=1 && size(x,2) ~= size(y,2)
        error(message('signal:computeperiodogram:MismatchedNumberOfChannels'))
    end
end

if Lx ~= Lw,
    error(message('signal:computeperiodogram:invalidWindow', 'WINDOW'))
end

if (numel(x)<2 || numel(size(x))>2)
    error(message('signal:computeperiodogram:NDMatrixUnsupported'))
end
