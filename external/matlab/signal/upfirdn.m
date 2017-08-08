function Y = upfirdn(x,h,varargin)
%UPFIRDN  Upsample, apply a specified FIR filter, and downsample a signal.
%   UPFIRDN(X,H,P,Q) is a cascade of three systems applied to input signal X:
%         (1) Upsampling by P (zero insertion).  P defaults to 1 if not 
%             specified.
%         (2) FIR filtering with the filter specified by the impulse response 
%             given in H.
%         (3) Downsampling by Q (throwing away samples).  Q defaults to 1 if not 
%             specified.
%   UPFIRDN uses an efficient polyphase implementation.
%
%   Usually X and H are vectors, and the output is a (signal) vector. 
%   UPFIRDN permits matrix arguments under the following rules:
%   If X is a matrix and H is a vector, each column of X is filtered through H.
%   If X is a vector and H is a matrix, each column of H is used to filter X.
%   If X and H are both matrices with the same number of columns, then the i-th
%      column of H is used to filter the i-th column of X.
%
%   Specifically, these rules are carried out as follows.  Note that the length
%   of the output is Ly = ceil( ((Lx-1)*P + Lh)/Q ) where Lx = length(X) and 
%   Lh = length(H). 
%
%      Input Signal X    Input Filter H    Output Signal Y   Notes
%      -----------------------------------------------------------------
%   1) length Lx vector  length Lh vector  length Ly vector  Usual case.
%   2) Lx-by-Nx matrix   length Lh vector  Ly-by-Nx matrix   Each column of X
%                                                            is filtered by H.
%   3) length Lx vector  Lh-by-Nh matrix   Ly-by-Nh matrix   Each column of H is
%                                                            used to filter X.
%   4) Lx-by-N matrix    Lh-by-N matrix    Ly-by-N matrix    i-th column of H is
%                                                            used to filter i-th
%                                                            column of X.
%
%   For an easy-to-use alternative to UPFIRDN, which does not require you to 
%   supply a filter or compensate for the signal delay introduced by filtering,
%   use RESAMPLE.
%
%   EXAMPLE: Sample-rate conversion by a factor of 147/160. It is used to
%            % downconvert from 48kHz to 44.1kHz.
%            L = 147; M = 160;                   % Interpolation/decimation factors.
%            Lp = 24;                            % Filter length of each phase
%            N = Lp*L-1;                         % Filter Order
%            h = fir1(N,1/M,kaiser(N+1,7.8562));
%            h = L*h; % Passband gain = L
%            Fs = 48e3;                          % Original sampling frequency: 48kHz
%            n = 0:10239;                        % 10240 samples, 0.213 seconds long
%            x  = sin(2*pi*1e3/Fs*n);            % Original signal, sinusoid at 1kHz
%            y = upfirdn(x,h,L,M);               % 9430 samples, still 0.213 seconds
%
%            % Overlay original (48kHz) with resampled signal (44.1kHz) in red.
%            stem(n(1:49)/Fs,x(1:49)); hold on 
%            stem(n(1:45)/(Fs*L/M),y(12:56),'r','filled'); 
%            xlabel('Time (sec)');ylabel('Signal value');
%    
%   See also RESAMPLE, INTERP, DECIMATE, FIR1, INTFILT.
  
%   Author(s): Paul Pacheco
%   Copyright 1988-2010 The MathWorks, Inc.

% This file validates the inputs, sets defaults, and then calls the C
% MEX-file.

% Validate number of I/O args.
error(nargchk(2,4,nargin,'struct'));
error(nargoutchk(0,1,nargout,'struct'));

% Force to be a column if input is a vector
[mx,nx] = size(x);
if find([mx nx]==1),
  x = x(:);  % columnize it.
end
[Lx,nChans] = size(x);

% Force to be a column if filter is a vector
if find(size(h)==1),
    h = h(:);  % columnize it.
end
[Lh,hCols] = size(h);

% Validate input args and define defaults.
[p,q] = validateinput(x,h,varargin);

% Call the MEX-file
Y = upfirdnmex(x,h,p,q,Lx,Lh,hCols,nChans);

% Convert output to be a row vector (if x was a row and H is NOT a matrix)
if (mx==1) && (hCols == 1)
    Y = Y(:).';
end


%----------------------------------------------------------------------
function [p,q] = validateinput(x,h,opts)

% Default values
p = 1;
q = 1;

% Validate 1st two input args: signal and filter.
if isempty(x) || issparse(x) || ~isa(x,'double'),
    error(message('signal:upfirdn:invalidInput', 'X'));
end
if isempty(h) || issparse(h) || ~isa(h,'double'),
    error(message('signal:upfirdn:invalidFilter', 'H'));
end

% The following check is for case 4 (as seen on the reference page), i.e., 
% x and h are matrices, check that they both have the same number of
% columns. 
nChans = size(x, 2);
hCols  = size(h, 2);
if (nChans > 1) && (hCols > 1) && (hCols ~= nChans),
    error(message('signal:upfirdn:xNhSizemismatch', 'X', 'H'));
end

% Validate optional input args: upsample and downsample factors.
nopts = length(opts);
if (nopts >= 1),
    p = opts{1};
    if isempty(p) || ~isa(p,'double') || p<1 || ~isequal(round(p),p),
        error(message('signal:upfirdn:invalidP', 'P'));
    elseif (nopts == 2),
        q = opts{2};
        if isempty(q) || ~isa(q,'double') || q<1 || ~isequal(round(q),q),
            error(message('signal:upfirdn:invalidQ', 'Q'));
        end
    end
    if p*q > intmax('int32'),
        error(message('signal:upfirdn:ProdPNQTooLarge', 'Q', 'P'));
    end
end

% [EOF]

